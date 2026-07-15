#!/usr/bin/env python3
"""Bridge benchtools assembly lists into an MGSIM genome table.

`makeBenchmarkSet` and `makeInclusionQuerySet` emit plain lists of assembly
accessions (one per line), e.g. the `.databaseAssembly`, `.totalExcludedAssembly`,
`.subspeciesInclusionAssemblies`, and `.speciesInclusionAssemblies` files. MGSIM
(https://github.com/nick-youngblut/MGSIM) simulates metagenomic samples from a
tab-separated genome table.

This script turns one or more accession lists into that table. MGSIM uses two
table flavours:

  * `MGSIM genome_download` expects `Taxon` + `Accession`.
  * `MGSIM communities` / `MGSIM reads` expect `Taxon` + `Fasta` (a path to a
    genome FASTA on disk).

Because the benchtools emit *assembly* accessions (GCF_/GCA_) while MGSIM's
`genome_download` expects *nucleotide* accessions (NC_/NZ_), the recommended path
is to point this script at the genome FASTAs you already have (the same files you
build the classifier database from) with --genome-dir or --fasta-map. It then
resolves each accession to its FASTA and writes a `Taxon` + `Accession` + `Fasta`
table that `communities`/`reads` can consume directly.

The `Taxon` label is the assembly accession itself, so every simulated read
traces straight back to its source assembly for grading against an
accession -> taxid answer sheet.

Examples
--------
Genome table pointing at local FASTAs (feeds communities/reads directly):

    ./benchtools_to_mgsim.py \
        assemblies.list.subspeciesInclusionAssemblies \
        --genome-dir /data/gtdb/genomes \
        -o mgsim.genomes.tsv

Accession-only table (feeds MGSIM genome_download, which downloads the FASTAs):

    ./benchtools_to_mgsim.py assemblies.list.totalExcludedAssembly \
        -o mgsim.taxon_accession.tsv

Then, with MGSIM installed (`pip install MGSIM`):

    MGSIM communities --n-comm 2 mgsim.genomes.tsv out/comm
    MGSIM reads mgsim.genomes.tsv --sr-seq-depth 1e6 out/comm_abund.txt out/reads/
"""

import argparse
import os
import sys
from glob import glob


def strip_version(accession):
    """GCF_000006945.2 -> GCF_000006945"""
    return accession.split(".", 1)[0]


def read_accessions(paths, column, has_header):
    """Yield (accession, source_path, lineno) from one-accession-per-line files.

    Each line is split on whitespace/comma and the chosen column is taken. Blank
    lines are skipped; a leading header line is skipped when --has-header is set.
    """
    for path in paths:
        with open(path) as handle:
            for lineno, raw in enumerate(handle, start=1):
                if has_header and lineno == 1:
                    continue
                line = raw.strip()
                if not line:
                    continue
                fields = line.replace(",", "\t").split("\t")
                if column >= len(fields):
                    sys.exit(
                        f"{path}:{lineno}: line has {len(fields)} column(s), "
                        f"cannot take column {column}"
                    )
                acc = fields[column].strip()
                if acc:
                    yield acc, path, lineno


def load_fasta_map(path):
    """Two-column TSV: accession <tab> fasta_path."""
    mapping = {}
    with open(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                sys.exit(f"{path}:{lineno}: expected 'accession<TAB>path'")
            mapping[parts[0].strip()] = parts[1].strip()
    return mapping


def resolve_from_dir(accession, genome_dir, suffixes):
    """Find a FASTA for `accession` under `genome_dir`.

    Tries the accession as given and with its version stripped, against each
    allowed suffix and a permissive glob. Returns the path or None.
    """
    candidates = [accession, strip_version(accession)]
    for acc in candidates:
        for suffix in suffixes:
            direct = os.path.join(genome_dir, acc + suffix)
            if os.path.isfile(direct):
                return direct
    for acc in candidates:
        hits = sorted(glob(os.path.join(genome_dir, acc + "*")))
        hits = [h for h in hits if os.path.isfile(h)]
        if hits:
            return hits[0]
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert benchtools assembly lists into an MGSIM genome table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "accession_lists",
        nargs="+",
        help="one or more benchtools output files (one accession per line)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="output genome table (default: stdout)",
    )
    parser.add_argument(
        "--column",
        type=int,
        default=0,
        help="0-based column holding the accession when lines have several "
        "fields (default: 0)",
    )
    parser.add_argument(
        "--has-header",
        action="store_true",
        help="skip the first line of each input file",
    )
    parser.add_argument(
        "--genome-dir",
        help="directory of genome FASTAs; each accession is resolved to a file "
        "here so the table feeds MGSIM communities/reads directly",
    )
    parser.add_argument(
        "--fasta-suffix",
        action="append",
        default=None,
        help="FASTA suffix(es) to try under --genome-dir (repeatable); "
        "default: several common .fna/.fasta[.gz] endings",
    )
    parser.add_argument(
        "--fasta-map",
        help="two-column TSV (accession<TAB>fasta_path) to resolve FASTAs "
        "instead of, or in addition to, --genome-dir",
    )
    parser.add_argument(
        "--no-dedup",
        action="store_true",
        help="keep duplicate accessions (default: drop duplicates, first wins)",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="warn instead of failing when a FASTA cannot be resolved; the row "
        "is written with an empty Fasta column",
    )
    args = parser.parse_args(argv)

    suffixes = args.fasta_suffix or [
        "_genomic.fna.gz",
        "_genomic.fna",
        ".fna.gz",
        ".fna",
        ".fasta.gz",
        ".fasta",
        ".fa.gz",
        ".fa",
    ]

    fasta_map = load_fasta_map(args.fasta_map) if args.fasta_map else {}
    want_fasta = bool(args.genome_dir or fasta_map)

    seen = set()
    rows = []
    missing = []
    for acc, path, lineno in read_accessions(
        args.accession_lists, args.column, args.has_header
    ):
        if not args.no_dedup:
            if acc in seen:
                continue
            seen.add(acc)

        fasta = None
        if want_fasta:
            fasta = fasta_map.get(acc) or fasta_map.get(strip_version(acc))
            if not fasta and args.genome_dir:
                fasta = resolve_from_dir(acc, args.genome_dir, suffixes)
            if not fasta:
                missing.append((acc, path, lineno))
        rows.append((acc, fasta))

    if missing and not args.allow_missing:
        for acc, path, lineno in missing[:20]:
            print(f"error: no FASTA found for {acc} ({path}:{lineno})", file=sys.stderr)
        if len(missing) > 20:
            print(f"error: ... and {len(missing) - 20} more", file=sys.stderr)
        sys.exit(
            f"error: {len(missing)} accession(s) could not be resolved to a FASTA. "
            "Fix --genome-dir/--fasta-map/--fasta-suffix, or pass --allow-missing."
        )

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        if want_fasta:
            out.write("Taxon\tAccession\tFasta\n")
            for acc, fasta in rows:
                out.write(f"{acc}\t{acc}\t{fasta or ''}\n")
        else:
            out.write("Taxon\tAccession\n")
            for acc, _ in rows:
                out.write(f"{acc}\t{acc}\n")
    finally:
        if out is not sys.stdout:
            out.close()

    dest = "stdout" if args.output == "-" else args.output
    summary = f"Wrote {len(rows)} genome(s) to {dest}"
    if want_fasta and missing and args.allow_missing:
        summary += f" ({len(missing)} without a resolved FASTA)"
    print(summary, file=sys.stderr)


if __name__ == "__main__":
    main()

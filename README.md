# metabuli-benchtools

Standalone benchmark & grading tools extracted from
[Metabuli](https://github.com/steineggerlab/Metabuli). This project bundles the
three tools that are useful on their own for building benchmark sets and scoring
taxonomic classification results, **without depending on the full mmseqs2
library** — so it compiles in seconds instead of minutes.

## Tools

| Command | Purpose |
| --- | --- |
| `grade` | Score classification results against answer sheets (precision / sensitivity / F1 per rank). |
| `makeBenchmarkSet` | Build family/genus/species/subspecies exclusion sets, the database assembly list, and species/subspecies inclusion query pairs from an assembly list (GTDB or virus). |
| `sample-queries` | Draw a diversity-maximizing subset of query genomes from a `makeBenchmarkSet` manifest, with per-category ratios. |

## Build

Requires a C++17 compiler and CMake ≥ 3.10. OpenMP is optional (enables
multi-threaded grading over multiple input files).

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The binary is produced at `build/benchtools`.

## Usage

```sh
benchtools <command> [options]
benchtools <command> --help     # command-specific usage
```

### grade

```sh
benchtools grade <classificationList> <mappingList> <taxonomyDir> [options]
```

- `<classificationList>` — file listing one classification-result path per line.
- `<mappingList>` — file listing one accession→taxid answer-sheet path per line
  (aligned line-by-line with the classification list).
- `<taxonomyDir>` — directory containing `names.dmp`, `nodes.dmp`, `merged.dmp`.

Key options: `--test-type` (`gtdb` [default], `gtdb-amgsim`, `cami`, `cami-long`,
`cami-euk`, `hiv`, `hiv-ex`, `over`, `kapk`), `--rank`, `--read-id-col`,
`--tax-id-col`, `--print-cols`, `--skip-secondary`, `--threads`.

### makeBenchmarkSet

```sh
benchtools makeBenchmarkSet <assemblyList> <taxonomyDir> [--test-type gtdb|virus] [--seed N] [--acc2taxid FILE] [--prefix STR] [--skip-validation]
```

`--acc2taxid` is required for `--test-type virus`. `--skip-validation` skips the
per-query exclusion/inclusion sanity checks (faster on large sets; the output
files are unchanged). `--prefix` sets the output path prefix (default:
`<assemblyList>`).

A single run writes three files, named `<prefix>.*` (with `<prefix>` defaulting to
`<assemblyList>`):

- **`.database`** — the assemblies to build the classifier database from
  (everything not held out), one accession per line.
- **`.summary`** — totals (input / database / excluded genome counts and the
  number of orders/families/genera/species in the input tree) plus a per-rank
  table of exclusion/inclusion counts.
- **`.query.tsv`** — one row per query assembly, with the answer-key metadata
  needed to grade it. Columns:

  | Column | Meaning |
  | --- | --- |
  | `Accession` | the query assembly |
  | `Category` | which test this query belongs to (see below) |
  | `ExpectedRank` | deepest rank a correct call can reach given the database |
  | `QueryTaxID` | taxid of the query assembly itself |
  | `SubjectTaxID` | the taxon defining the case (excluded/shared); paired rows share it |
  | `SubjectRank` | rank of `SubjectTaxID` |
  | `SubjectSize` | number of direct members of `SubjectTaxID` |

  `Category` is one of `familyExclusion`, `genusExclusion`, `speciesExclusion`,
  `subspeciesExclusion` (a held-out taxon whose siblings remain, so a correct
  call reaches `ExpectedRank` = order/family/genus/species respectively); or
  `speciesInclusionPair` / `subspeciesInclusionPair` (GTDB only — two species of
  the same genus / two genomes of the same species, sharing `SubjectTaxID`).

For `familyExclusion`, `genusExclusion`, and `speciesExclusion`, the manifest
lists **every member genome** of each excluded taxon (the whole taxon is out of
the database, so all its genomes are valid queries). Rows are grouped by
`SubjectTaxID` so a subset can be sampled afterward with `sample-queries` (below).

> Earlier versions wrote ~11 files and had a separate `makeInclusionQuerySet`
> command; both are folded into this single command and its `.database` /
> `.query.tsv` / `.summary` output. The
> database composition is unchanged; the exclusion query sets now enumerate all
> excluded genomes (earlier versions reported one sampled query per excluded
> taxon), and the single per-species inclusion set was dropped in favor of the
> inclusion pairs.

### sample-queries

```sh
benchtools sample-queries <queryTsv> <outPrefix> --number N [--ratio a,b,c,d,e,f] [--seed S]
```

Draws a subset of query genomes from a `makeBenchmarkSet` `.query.tsv` manifest.
`<outPrefix>` is required; the tool writes two files:

- **`<outPrefix>.query.tsv`** — the sampled rows, same columns as the input.
- **`<outPrefix>.summary`** — the input path, the sampling parameters
  (`number` / `ratio` / `seed`), a per-category sampling table
  (`Available` / `Requested` / `Sampled`), and any capping warnings.

`--number` is the total sample size and `--ratio` is six comma-separated weights in
the order `familyExcl,genusExcl,speciesExcl,subspeciesExcl,speciesIncl,subspeciesIncl`
(default all `1`) that split the total across categories.

Within each category the sampler **maximizes diversity across taxa** by
round-robin over `SubjectTaxID` groups — because `makeBenchmarkSet` excludes one
taxon per parent, this spreads picks over as many distinct parents as possible
(e.g. species-exclusion queries come from as many different genera as available).

The two inclusion categories are pairs: their weight counts **pairs**, and both
members of a chosen pair are always kept together (each pair adds two rows). If a
category has fewer units than its target, all are taken and a warning is recorded
(in the summary and on stderr), so the total may fall short of `--number`.

```sh
benchtools sample-queries assemblies.list.query.tsv sample1 --number 300 --ratio 1,2,3,4,5,6 --seed 0
# -> sample1.query.tsv (feed its column 0 to the MGSIM bridge below) + sample1.summary
```

## Simulating metagenomic samples with MGSIM

`makeBenchmarkSet` writes `<assemblyList>.query.tsv`, whose first column is the
assembly accession of every genome from which benchmark reads should be
simulated (filter by `Category` for a specific test). `<assemblyList>.database`
lists the accessions that instead go into the classifier database.

[MGSIM](https://github.com/nick-youngblut/MGSIM) turns a genome table into
simulated Illumina / PacBio / Nanopore metagenomes. It is an external Python
tool — install it separately:

```sh
pip install MGSIM        # or: conda env create -f MGSIM/environment.yml
```

`scripts/benchtools_to_mgsim.py` bridges a benchtools accession list into the
tab-separated genome table MGSIM expects. It uses the **assembly accession as the
`Taxon` label**, so every simulated read traces straight back to its source
assembly for grading against an accession→taxid answer sheet (`benchtools grade`).

MGSIM uses two table flavours, and the bridge produces whichever you need:

- **`Taxon` + `Fasta`** (for `MGSIM communities` / `MGSIM reads`). Point the
  bridge at the genome FASTAs you already have — usually the same files used to
  build the classifier database — and it resolves each accession to a file:

  ```sh
  scripts/benchtools_to_mgsim.py \
      assemblies.list.query.tsv --has-header \
      --genome-dir /data/gtdb/genomes \
      -o mgsim.genomes.tsv

  MGSIM communities --n-comm 2 mgsim.genomes.tsv out/comm
  MGSIM reads mgsim.genomes.tsv --sr-seq-depth 1e6 out/comm_abund.txt out/reads/
  ```

  (The bridge reads column 0 by default; `--has-header` skips the manifest's
  header row. To simulate only one test, filter first, e.g.
  `awk -F'\t' '$2=="familyExclusion"{print $1}' assemblies.list.query.tsv`.)

- **`Taxon` + `Accession`** (for `MGSIM genome_download`, which fetches the
  FASTAs itself). Omit `--genome-dir`:

  ```sh
  scripts/benchtools_to_mgsim.py assemblies.list.query.tsv --has-header \
      -o mgsim.taxon_accession.tsv
  MGSIM genome_download -d genomes/ mgsim.taxon_accession.tsv > mgsim.genomes.tsv
  ```

  Note: benchtools emit *assembly* accessions (`GCF_`/`GCA_`) while MGSIM's
  `genome_download` expects *nucleotide* accessions (`NC_`/`NZ_`). If you already
  have the assembly FASTAs, prefer the `--genome-dir` path above.

Run `scripts/benchtools_to_mgsim.py --help` for all options (`--fasta-map`,
`--fasta-suffix`, `--column`, `--has-header`, `--allow-missing`, `--no-dedup`).

## Relationship to Metabuli

The taxonomy engine (`NcbiTaxonomy`, `TaxonomyWrapper`) and the tool
sources are ported from Metabuli. Everything the tools needed from the mmseqs2
framework has been replaced by small self-contained shims under `src/compat/`:

- **CLI** — Metabuli's `LocalParameters` / mmseqs `Parameters` / `Command`
  framework is replaced by a lightweight `Parameters` struct and hand-rolled
  argument parser (`src/Parameters.*`, `src/main.cpp`).
- **Utilities** — `Util::split`, `FileUtil::fileExists`, `MathUtil::flog2/ipow`,
  `Debug`/`EXIT`, and `SORT_SERIAL` are reimplemented in `src/compat/`.
- Unused dependencies (`IndexCreator`, `KSeqWrapper`, `DBReader`, …) were dropped.

The tools produce output byte-identical to the corresponding Metabuli expert
commands (`grade`, and `maketestsets` + `makeInclusionTestQueries` combined —
the latter two are folded into the single `makeBenchmarkSet` command here).

## Layout

```
src/
  main.cpp            CLI dispatcher
  Parameters.{h,cpp}  lightweight options + argument parser
  Assembly.h          the Assembly struct (from Metabuli's common.h)
  compat/             minimal drop-in replacements for mmseqs headers
  taxonomy/           NcbiTaxonomy + TaxonomyWrapper (ported)
  tools/              grade, makeBenchmarkSet
scripts/
  benchtools_to_mgsim.py   accession list -> MGSIM genome table
```

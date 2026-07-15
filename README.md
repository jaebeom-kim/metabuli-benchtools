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
| `makeBenchmarkSet` | Build family/genus/species/subspecies exclusion + inclusion benchmark sets from an assembly list (GTDB or virus). |
| `makeInclusionQuerySet` | Build subspecies/species inclusion query sets from an assembly list. |

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
benchtools makeBenchmarkSet <assemblyList> <taxonomyDir> [--test-type gtdb|virus] [--seed N] [--acc2taxid FILE]
```

`--acc2taxid` is required for `--test-type virus`.

### makeInclusionQuerySet

```sh
benchtools makeInclusionQuerySet <assemblyList> <taxonomyDir>
```

## Relationship to Metabuli

The taxonomy engine (`NcbiTaxonomy`, `TaxonomyWrapper`) and the three tool
sources are ported from Metabuli. Everything the tools needed from the mmseqs2
framework has been replaced by small self-contained shims under `src/compat/`:

- **CLI** — Metabuli's `LocalParameters` / mmseqs `Parameters` / `Command`
  framework is replaced by a lightweight `Parameters` struct and hand-rolled
  argument parser (`src/Parameters.*`, `src/main.cpp`).
- **Utilities** — `Util::split`, `FileUtil::fileExists`, `MathUtil::flog2/ipow`,
  `Debug`/`EXIT`, and `SORT_SERIAL` are reimplemented in `src/compat/`.
- Unused dependencies (`IndexCreator`, `KSeqWrapper`, `DBReader`, …) were dropped.

The three tools produce output byte-identical to the corresponding Metabuli
expert commands (`grade`, `maketestsets`, `makeInclusionTestQueries`).

## Layout

```
src/
  main.cpp            CLI dispatcher
  Parameters.{h,cpp}  lightweight options + argument parser
  Assembly.h          the Assembly struct (from Metabuli's common.h)
  compat/             minimal drop-in replacements for mmseqs headers
  taxonomy/           NcbiTaxonomy + TaxonomyWrapper (ported)
  tools/              grade, makeBenchmarkSet, makeInclusionQuerySet
```

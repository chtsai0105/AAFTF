# CLAUDE.md — Claude Code Project Context for AAFTF

See [AGENTS.md](AGENTS.md) for full development guidelines, code style, and common patterns.

## Key Architecture Reminders

- Entry point: `AAFTF/AAFTF_main.py` — builds the top-level parser, wires up subcommands via `SUBCOMMAND_REGISTRARS` (imported from `AAFTF/_menu.py`), and dispatches by calling `args.func(**vars(args))`. There is no central dispatcher function: each subcommand parser binds its own subtool's `run` directly via `parser_x.set_defaults(func=<module>.run)` inside its `<name>_menu()` function in `_menu.py`. Subcommand parsers no longer register any `aliases=[...]` — each subcommand has exactly one canonical name.
- All argparse subcommand-parser definitions ("menus") live in one place: `AAFTF/_menu.py`. Each subcommand has a `<name>_menu(subparsers)` function there (e.g. `trim_menu`, `depth_menu`) that builds and registers that subtool's `subparsers.add_parser(...)` block and its arguments, then calls `parser_x.set_defaults(func=<module>.run)`; `AAFTF/_menu.py` also defines `menu_common_args()` and the `SUBCOMMAND_REGISTRARS` list, and imports every subcommand module (`trim`, `depth`, etc.) to reference their `run` functions. Each subcommand module itself (`trim.py`, `depth.py`, etc.) only contains the `run(**kwargs)` execution logic — no parser-building code
- Every subcommand's `run()` function takes its CLI dest names directly as keyword arguments (not an `args`/`parser` pair) and ends its parameter list with `**kwargs` to absorb argparse's own bookkeeping attributes (`command`, `func`, and `quiet`, which only `main()` uses to set the log level) that ride along in `vars(args)`. The function body uses those parameters directly (e.g. `left`, `workdir`) — there is no internal `args = Namespace(...)` reconstruction; a parameter that needs to change during the run (e.g. `basename` auto-derived from `left`, `workdir` defaulted when not given) is simply reassigned as a local variable. Only `assemble.py`'s internal helpers (`run_spades`, `run_dipspades`, `run_megahit`, `run_unicycler`) and `pipeline.py`'s `create_namespace()` still build/pass `Namespace`/dict objects, because they need a keyed collection to forward a subset of fields — not because they mirror the old `args.xxx` style.
- Because `run()` no longer mutates a caller-supplied `args` object, tests must not expect derived values (e.g. `trim.run()`'s auto-derived `basename`) to show up on the `Namespace`/kwargs the test constructed — assert on the actual side effects (subprocess commands, files written) instead. See `tests/test_trim.py::_run_bbduk` for the pattern.
- `AAFTF/pipeline.py`'s own `run()` follows the same convention, and internally calls each step's module via `<module>.run(**vars(step_args))` instead of `<module>.run(parser, step_args)`
- `CustomHelpFormatter` lives in `AAFTF/utility.py`. Every subcommand parser (built in `_menu.py`) uses `formatter_class=CustomHelpFormatter`

## CLI Framework Conventions

- Every subcommand parser must have `-v/--verbose` (argparse dest `debug`, so code and `run()` parameters still use `debug`) and `--pipe` flags — enforced by test suite (`TestAssessParser`, `TestSortParser`, `TestFixTblParser`)
- When building a `Namespace` for pipeline steps, include **all** attributes the target `run()` function accesses — check the submodule source to avoid `AttributeError`
- `menu_common_args(target)` (in `AAFTF/_menu.py`) adds exactly three arguments — `--pipe`, `-q/--quiet`, `-v/--verbose` (dest `debug`) — to whatever `target` is passed. Every `<name>_menu()` function calls it **last**, passing its own `optional = parser_x.add_argument_group("optional arguments")` group (not the parser itself), so these three common flags render as the final entries of that subcommand's "optional arguments" section rather than in a separate leading group. Do not pass the raw parser to `menu_common_args()` — always pass the `optional` group, and call it after all of that subcommand's own optional args have been added.
- All other args a subcommand needs (`-c/--cpus`, `-w/--workdir/--tmpdir`, `-l/--left`/`-r/--right`, etc.) are implemented locally within that subcommand's own `<name>_menu()` in `_menu.py` — there is no shared parent-parser mechanism for them (an earlier version used `argparse` `parents=[...]` for these, but that shared the underlying `Action` objects across every subcommand and any per-subcommand override via `conflict_handler="resolve"` silently corrupted other subcommands; it was removed for this reason). When adding a `-c/--cpus`/`-w/--workdir`/`-l/--left`,`-r/--right` block to a new menu function, copy it verbatim from a similar existing one (e.g. `vecscreen_menu`, `sourpurge_menu`) to keep help text/defaults consistent.
- Each `<name>_menu()` function that has any required argument also creates two argument groups — `required = parser_x.add_argument_group("required arguments")` and `optional = parser_x.add_argument_group("optional arguments")` — and adds each local argument to the appropriate one (required = `required=True` with no `default=`; everything else = optional).
- Adding a new subcommand: write its `<name>_menu(subparsers)` function in `AAFTF/_menu.py` (create the subparser, add required/optional groups and any subtool-specific args, call `menu_common_args(optional)` last, then `parser_x.set_defaults(func=<module>.run)`), add it to `SUBCOMMAND_REGISTRARS` in that same file, and write the subtool's `run(**kwargs)` in its own module (with named params for the dest names it uses, plus a trailing `**kwargs`)
- Testing a subcommand's CLI parsing without executing the tool: patch `AAFTF.<module>.run` (not a central dispatcher — there isn't one) with a `side_effect` that captures `**kwargs` into a `Namespace`, then call `AAFTF.AAFTF_main.main()` with `sys.argv` patched. See `tests/test_cli.py::_parse_with_main` for the pattern. Note that since `run()` rebuilds its own internal `args`/local variables from the passed kwargs, a test's original `Namespace`/kwargs object is **not** mutated by any in-`run()` auto-derivation (e.g. `trim.run()` deriving `basename` from `left` when not given) — assert on side effects (subprocess commands, files written) instead of re-inspecting the caller's original args object for such derived values.

## Subcommand Reference

"Module" below is where each subtool's `run(**kwargs)` execution logic lives; its argparse parser/menu is instead in `AAFTF/_menu.py` (as `<name>_menu()`).

| Canonical name | Module | Key external tools |
|---|---|---|
| `trim` | `AAFTF/trim.py` | bbduk.sh, trimmomatic, fastp |
| `filter` | `AAFTF/filter.py` | bbduk.sh, bowtie2, bwa, minimap2, samtools |
| `assemble` | `AAFTF/assemble.py` | spades.py, megahit, unicycler |
| `vecscreen` | `AAFTF/vecscreen.py` | blastn, makeblastdb |
| `fcs_screen` | `AAFTF/fcs_screen.py` | run_fcsadaptor.sh (singularity/docker) |
| `fcs_gx_purge` | `AAFTF/fcs_gx_purge.py` | run_gx.py (NCBI FCS-GX) |
| `sourpurge` | `AAFTF/sourpurge.py` | sourmash, bwa, samtools |
| `rmdup` | `AAFTF/rmdup.py` | minimap2 |
| `polish` | `AAFTF/polish.py` | polypolish, pypolca, nextPolish2, racon, bwa, samtools, minimap2, yak, freebayes |
| `sort` | `AAFTF/sort.py` | (none — BioPython only) |
| `assess` | `AAFTF/assess.py` | (none — BioPython only) |
| `depth` | `AAFTF/depth.py` | minimap2, bwa, samtools, mosdepth |
| `mito` | `AAFTF/mito.py` | NOVOPlasty, minimap2 |
| `fix_tbl` | `AAFTF/fix_tbl.py` | (none) |
| `download` | `AAFTF/download.py` | (none — urllib only) |
| `check_dependencies` | `AAFTF/check_dependencies.py` | (none — checks PATH for all of the above) |
| `pipeline` | `AAFTF/pipeline.py` | all of the above |

## Module `args` Signatures

Key `args` attributes accessed by each `run()` function:

- **trim**: `left`, `right`, `basename`, `method`, `memory`, `cpus`, `minlen`, `avgqual`, `debug`, `pipe`, `trimmomatic`, `trimmomatic_adaptors`, `trimmomatic_leadingwindow`, `trimmomatic_trailingwindow`, `trimmomatic_slidingwindow`, `trimmomatic_quality`, `trimmomatic_clip`, `merge`, `dedup`, `cutfront`, `cuttail`, `cutright`
- **filter**: `workdir`, `cpus`, `left`, `right`, `screen_accessions`, `screen_urls`, `screen_local`, `basename`, `aligner`, `memory`, `debug`, `pipe`
- **assemble**: `method`, `workdir`, `cpus`, `memory`, `isolate`, `careful`, `assembler_args`, `tmpdir`, `left`, `right`, `merged`, `pipe`, `debug`, `out`
- **vecscreen**: `workdir`, `infile`, `outfile`, `cpus`, `percent_id`, `stringency`, `debug`, `pipe`
- **fcs_screen**: `container_engine`, `workdir`, `infile`, `image`, `prok`, `fcs_script`, `outfile`, `debug`, `pipe`
- **fcs_gx_purge**: `workdir`, `db`, `input`, `taxid`, `outfile`, `debug`, `pipe`
- **sourpurge**: `workdir`, `cpus`, `left`, `right`, `sourdb`, `sourdb_type`, `input`, `kmer`, `phylum`, `mincovpct`, `outfile`, `taxonomy`, `debug`, `pipe`
- **rmdup**: `workdir`, `cpus`, `input`, `percent_id`, `percent_cov`, `minlen`, `exhaustive`, `debug`, `out`, `pipe`
- **polish**: `method`, `memory`, `cpus`, `left`, `right`, `longreads`, `workdir`, `infile`, `outfile`, `debug`, `pipe`. Each `--method` (`pypolca`, `polypolish`, `nextpolish2`, `racon`) is isolated into its own `run_<method>()` function in `AAFTF/polish.py`.
- **sort**: `input`, `minlen`, `out`, `name`
- **assess**: `input`, `report`, `telomere_monomer`, `telomere_n_repeat`, `telomere_window`
- **depth**: `input`, `out`, `left`, `right`, `longreads`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `pipe`, `min_contig_len`, `no_plot`, `plot_format`
- **mito**: `workdir`, `left`, `right`, `seed`, `reference`, `minlen`, `maxlen`, `out`, `starting`, `memory`, `debug`, `pipe`
- **fix_tbl**: `table`, `report`, `output`, `debug`, `pipe`

## Pipeline `create_namespace()` Convention

`pipeline.py` uses a `create_namespace(options, required_args=None, **extra)` helper:
- `options`: list of keys to copy from the pipeline `args_dict`
- `required_args`: dict of values to override/add (pipeline-step-specific settings)
- Always include `"pipe": True` in `required_args` for every step
- Always include `"debug"` in the `options` list so parent debug flag propagates
- Do NOT append attributes directly to the returned Namespace; include them in `required_args` before calling

## Pipeline Workflow

Recommended order for full assembly QC:

```
trim → [mito optional, PE only] → filter → assemble → vecscreen
  → [fcs_screen / fcs_gx_purge optional] → sourpurge → rmdup
  → polish → sort → assess → [depth optional]
```

Expected inputs/outputs per step:

| Step | Input | Output |
|---|---|---|
| trim | raw FASTQ | `{base}_1P.fastq.gz`, `{base}_2P.fastq.gz` |
| filter | trimmed FASTQ | `{base}_filtered_1.fastq.gz`, `{base}_filtered_2.fastq.gz` |
| assemble | filtered FASTQ | `{base}.{method}.fasta` |
| vecscreen | assembled FASTA | `{base}.vecscreen.fasta` |
| sourpurge | vecscreen FASTA | `{base}.sourpurge.fasta` |
| rmdup | sourpurge FASTA | `{base}.rmdup.fasta` |
| polish | rmdup FASTA | `{base}.polish.fasta` |
| sort | polished FASTA | `{base}.final.fasta` |
| assess | sorted FASTA | printed stats (optional report file) |
| depth | final FASTA + reads | `coverage_stats.txt` |

## External Tool Versions (minimum known-good)

- samtools >= 1.13 (pixi pins >= 1.24; needed for `sort --write-index` and `flagstat -O tsv`). No samtools version branching: aligner → sorted BAM pipelines use `utility.align_to_sorted_bam()`, which pipes SAM straight into `samtools sort` and writes the `.bai` index in the same step
- minimap2 >= 2.17
- mosdepth >= 0.3
- sourmash >= 4.x (LCA-based classification)
- SPAdes >= 3.15 for `assemble`
- NCBI FCS-adaptor 0.5.5 (hard-coded in `resources.py`)

## Recent Additions (2026-03-02)

### `depth` subtool — coverage analysis

**New file:** `AAFTF/depth.py`; registered in `AAFTF_main.py`.

**Workflow:**
1. Counts reads in each input FASTQ
2. Maps Illumina reads with `minimap2 -ax sr` (default) or `bwa mem` (`--aligner bwa`); maps long reads with `minimap2 -ax map-ont` (default) / `map-pb` / `map-hifi` (`--longread_type`)
3. Merges BAMs when both read types are provided (`samtools merge`)
4. Runs `samtools flagstat` per read type, then `mosdepth` on the combined BAM
5. Parses `mosdepth.summary.txt` for per-contig mean depths and `mosdepth.global.dist.txt` for coverage breadth
6. Flags contigs with mean depth > assembly_mean + 3×SD as possible contaminants/organelles (uses population SD; contigs are the full assembly population, not a sample)

**Key `args` attributes:** `input`, `out`, `left`, `right`, `longreads`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `pipe`, `min_contig_len`, `no_plot`, `plot_format`

## Audit Fixes (2026-05-02)

### Bugs fixed

| File | Line | Issue | Fix |
|---|---|---|---|
| `filter.py` | 225, 251, 273 | `tempfiles[3]` IndexError — only `tempfiles[0]` existed | Replaced with `unsorted_bam` path variable |
| `polish.py` | 82 | `polish_log = None` inside loop before `open()` | Removed erroneous reset; branches already set it |
| `polish.py` | 71 | `i` undefined if `--iterations 0` | Added `sys.exit(1)` guard when `iterations < 1` |
| `trim.py` | 135 | Infinite loop: `dirname("/") == "/"` | Added `if new_path == findpath: break` guard |
| `vecscreen.py` | 288 | Stale `start`/`end` in else-branch for multi-hit contigs | Moved `start, end = sorted(...)` before the if/else |
| `vecscreen.py` | 262 | `global contigs_to_remove` leaked state across calls | Changed to local variable in `run()` |
| `sourpurge.py` | 113 | `"nomatch" in cols` tested list membership — fragile | Changed to `cols[1].strip() == "nomatch"` |
| `fcs_screen.py` | 31 | Missing DB fell back to string literal `"AAFTF_DB"` | Changed to `status()` + `sys.exit(1)` |

### CLI consistency fixes (AAFTF_main.py)

- Added `--pipe` and `-v/--debug` to `fix_tbl` parser
- Added `--pipe` to `pipeline` parser
- Removed dead `runall` stub from `run_subtool()`
- Fixed assembler help text: "Default: unicycler" → "Default: spades"

### pipeline.py refactor

- `sort_args` and `assess_args` now use `create_namespace()` (consistent with all other steps)
- Assembly-specific attributes (`merged`, `isolate`, `careful`) now passed via `required_args` before `create_namespace()` call

### Bioinformatics enhancements

- `assess.py`: `findTelomere()` window is now configurable (`--telomere_window`, default 200 bp); added N-gap count and total N bases to stats output
- `rmdup.py`: now computes and reports both N50 and N75; filtering threshold uses N75 (as before)
- `depth.py`: added `--min_contig_len` (default 500 bp) to exclude short contigs from outlier depth analysis

## Running Tests

```bash
# Unit tests only (no external tools required)
conda run -n base python -m pytest tests/ -m unit -v

# All tests (integration tests require external tool environment)
conda run -n base python -m pytest tests/ -v

# Compile-check all modules
python -m py_compile AAFTF/filter.py AAFTF/polish.py AAFTF/trim.py \
    AAFTF/vecscreen.py AAFTF/sourpurge.py AAFTF/fcs_screen.py \
    AAFTF/depth.py AAFTF/AAFTF_main.py AAFTF/pipeline.py \
    AAFTF/assess.py AAFTF/rmdup.py AAFTF/_menu.py
```

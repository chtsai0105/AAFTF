# CLAUDE.md — Claude Code Project Context for AAFTF

See [AGENTS.md](AGENTS.md) for full development guidelines, code style, and common patterns.

## Key Architecture Reminders

- Entry point: `aaftf/main.py` — builds the top-level parser, wires up subcommands via `register_subcommands(parser)` (imported from `aaftf/_menu.py`; it creates the subparsers action and lists the subcommands in `AAFTF --help` under three group titles — Setup, Assembly pipeline, Annotation — using the help-only `SubcommandGroup` action from `aaftf/utility.py`), and dispatches by calling `args.func(**vars(args))`; `main()` returns a numeric exit code (0 ok, 1 error, 2 `FileNotFoundError`, 130 Ctrl-C) and subcommands raise exceptions instead of calling `sys.exit` (see AGENTS.md). There is no central dispatcher function: each subcommand parser binds its own subtool's `run` directly via `parser_x.set_defaults(func=<module>.run)` inside its `<name>_menu()` function in `_menu.py`. Subcommand parsers no longer register any `aliases=[...]` — each subcommand has exactly one canonical name.
- All argparse subcommand-parser definitions ("menus") live in one place: `aaftf/_menu.py`. Each subcommand has a `<name>_menu(subparsers)` function there (e.g. `trim_menu`, `depth_menu`) that builds and registers that subtool's `subparsers.add_parser(...)` block and its arguments, then calls `parser_x.set_defaults(func=<module>.run)`; `aaftf/_menu.py` also defines `add_verbosity_args()` and `register_subcommands()` (which calls every `<name>_menu()` group by group, in `AAFTF --help` order), and imports every subcommand module (`trim`, `depth`, etc.) to reference their `run` functions. Each subcommand module itself (`trim.py`, `depth.py`, etc.) only contains the `run(**kwargs)` execution logic — no parser-building code
- Every subcommand's `run()` function takes its CLI dest names directly as keyword arguments (not an `args`/`parser` pair) and ends its parameter list with `**kwargs` to absorb argparse's own bookkeeping attributes (`command`, `func`, and `quiet`, which only `main()` uses to set the log level) that ride along in `vars(args)`. The function body uses those parameters directly (e.g. `read1`, `workdir`) — there is no internal `args = Namespace(...)` reconstruction; a parameter that needs to change during the run (e.g. `basename` auto-derived from `read1`, `workdir` defaulted when not given) is simply reassigned as a local variable. Only `assemble.py`'s internal helpers (`run_spades`, `run_megahit`, `run_unicycler`) still build/pass `Namespace` objects, because they need a keyed collection to forward a subset of fields — not because they mirror the old `args.xxx` style.
- Because `run()` no longer mutates a caller-supplied `args` object, tests must not expect derived values (e.g. `trim.run()`'s auto-derived `basename`) to show up on the `Namespace`/kwargs the test constructed — assert on the actual side effects (subprocess commands, files written) instead. See `tests/test_trim.py::_run_bbduk` for the pattern.
- `aaftf/pipeline.py`'s own `run()` follows the same convention, and internally calls each step's module via `<module>.run(**vars(step_args))` instead of `<module>.run(parser, step_args)`
- `CustomHelpFormatter` lives in `aaftf/utility.py`. Every subcommand parser (built in `_menu.py`) uses `formatter_class=CustomHelpFormatter`

## CLI Framework Conventions

- Every subcommand parser must have `-v/--verbose` (argparse dest `debug`, so code and `run()` parameters still use `debug`) and `-q/--quiet` flags — enforced by test suite (`TestAssessParser`, `TestSortParser`, `TestFixTblParser`)
- A step that suggests the next command logs it with `logger.info("Your next command might be: ...")` inside `if not pipe:`. `pipe: bool = False` is a `run()`-only parameter (there is no `--pipe` CLI option — the one exception to "`run()` parameters are the CLI dest names"): the CLI never sets it, so hints show unless `-q/--quiet`; `pipeline._step_kwargs()` sets `pipe=True` so pipeline runs skip them.
- Log files: `setup_logging()` (called by `main()`) holds log messages in memory; `make_workdir(workdir, prefix)` writes them to `<workdir>/<prefix>.log` (append mode; `prefix` is the subcommand name) and keeps logging there, and `cleanup_workdir()` closes it before deciding whether to delete the workdir — so the log is kept exactly when the workdir is (own `-w`, `-v`, or a failed run). `main()` calls `finish_logging()` at the end: leftover messages (e.g. the next-command hint) are appended to that log if it still exists, and subcommands without a workdir write `./<command>.log` only with `-v`. The file gets INFO and above even with `-q` (DEBUG with `-v`) and always includes tracebacks; errors are logged once with `exc_info=True`, and the terminal shows the traceback only with `-v`. Always create workdirs with `make_workdir`/`cleanup_workdir` so steps get a log.
- `add_verbosity_args(target)` (in `aaftf/_menu.py`) adds exactly two arguments — `-q/--quiet` and `-v/--verbose` (dest `debug`) — to whatever `target` is passed. Every `<name>_menu()` function calls it **last**, passing its own `optional = parser_x.add_argument_group("optional arguments")` group (not the parser itself), so these two common flags render as the final entries of that subcommand's "optional arguments" section rather than in a separate leading group. Do not pass the raw parser to `add_verbosity_args()` — always pass the `optional` group, and call it after all of that subcommand's own optional args have been added.
- All other args a subcommand needs (`-c/--cpus`, `-w/--workdir/--tmpdir`, `-1/--read1`/`-2/--read2`, etc.) are implemented locally within that subcommand's own `<name>_menu()` in `_menu.py` — there is no shared parent-parser mechanism for them (an earlier version used `argparse` `parents=[...]` for these, but that shared the underlying `Action` objects across every subcommand and any per-subcommand override via `conflict_handler="resolve"` silently corrupted other subcommands; it was removed for this reason). When adding a `-c/--cpus`/`-w/--workdir`/`-1/--read1`,`-2/--read2` block to a new menu function, copy it verbatim from a similar existing one (e.g. `vecscreen_menu`, `sourpurge_menu`) to keep help text/defaults consistent.
- Path-like arguments (inputs, outputs, reads, databases, work dirs, basenames/prefixes) use `type=str`, never `argparse.FileType` or `pathlib.Path`: `run()` receives plain path strings and opens files itself (with `with open(...)`), using `Path` internally only where path operations help. `FileType` opens (and for outputs truncates) files at parse time and never closes them, which is why `fix_tbl` moved off it.
- Each `<name>_menu()` function that has any required argument also creates two argument groups — `required = parser_x.add_argument_group("required arguments")` and `optional = parser_x.add_argument_group("optional arguments")` — and adds each local argument to the appropriate one (required = `required=True` with no `default=`; everything else = optional).
- Adding a new subcommand: write its `<name>_menu(subparsers)` function in `aaftf/_menu.py` (create the subparser, add required/optional groups and any subtool-specific args, call `add_verbosity_args(optional)` last, then `parser_x.set_defaults(func=<module>.run)`), add it to the right group's list in `register_subcommands()` in that same file, and write the subtool's `run(**kwargs)` in its own module (with named params for the dest names it uses, plus a trailing `**kwargs`)
- Testing a subcommand's CLI parsing without executing the tool: patch `aaftf.<module>.run` (not a central dispatcher — there isn't one) with a `side_effect` that captures `**kwargs` into a `Namespace`, then call `aaftf.main.main()` with `sys.argv` patched. See `tests/test_cli.py::_parse_with_main` for the pattern. Note that since `run()` rebuilds its own internal `args`/local variables from the passed kwargs, a test's original `Namespace`/kwargs object is **not** mutated by any in-`run()` auto-derivation (e.g. `trim.run()` deriving `basename` from `read1` when not given) — assert on side effects (subprocess commands, files written) instead of re-inspecting the caller's original args object for such derived values.

## Subcommand Reference

"Module" below is where each subtool's `run(**kwargs)` execution logic lives; its argparse parser/menu is instead in `aaftf/_menu.py` (as `<name>_menu()`).

| Canonical name | Module | Key external tools |
|---|---|---|
| `trim` | `aaftf/trim.py` | bbduk.sh, trimmomatic, fastp |
| `filter` | `aaftf/filter.py` | bbduk.sh, bowtie2, bwa, minimap2, samtools |
| `assemble` | `aaftf/assemble.py` | spades.py, megahit, unicycler |
| `vecscreen` | `aaftf/vecscreen.py` | blastn, makeblastdb |
| `fcs_screen` | `aaftf/fcs_screen.py` | run_fcsadaptor.sh (singularity/docker) |
| `fcs_gx_purge` | `aaftf/fcs_gx_purge.py` | run_gx.py (NCBI FCS-GX) |
| `sourpurge` | `aaftf/sourpurge.py` | sourmash, bwa, samtools |
| `rmdup` | `aaftf/rmdup.py` | minimap2 |
| `polish` | `aaftf/polish.py` | polypolish, pypolca, nextPolish2, racon, bwa, samtools, minimap2, yak, freebayes |
| `sort` | `aaftf/sort.py` | (none — BioPython only) |
| `assess` | `aaftf/assess.py` | (none — BioPython only) |
| `depth` | `aaftf/depth.py` | minimap2, bwa, samtools, mosdepth |
| `mito` | `aaftf/mito.py` | NOVOPlasty, minimap2 |
| `fix_tbl` | `aaftf/fix_tbl.py` | (none) |
| `database` | `aaftf/database.py` | (none — urllib only) |
| `dependency` | `aaftf/dependency.py` | (none — checks PATH for all of the above) |
| `pipeline` | `aaftf/pipeline.py` | all of the above |

## Module `args` Signatures

Key `args` attributes accessed by each `run()` function:

- **trim**: `read1`, `read2`, `basename`, `method`, `memory`, `cpus`, `minlen`, `avgqual`, `debug`, `trimmomatic`, `trimmomatic_adaptors`, `trimmomatic_leadingwindow`, `trimmomatic_trailingwindow`, `trimmomatic_slidingwindow`, `trimmomatic_quality`, `trimmomatic_clip`, `merge`, `dedup`, `cutfront`, `cuttail`, `cutright`
- **filter**: `workdir`, `cpus`, `read1`, `read2`, `screen_accessions`, `screen_urls`, `screen_local`, `basename`, `aligner`, `memory`, `debug`
- **assemble**: `method`, `workdir`, `cpus`, `memory`, `isolate`, `careful`, `assembler_args`, `tmpdir`, `read1`, `read2`, `merged`, `debug`, `out`
- **vecscreen**: `workdir`, `infile`, `outfile`, `cpus`, `percent_id`, `stringency`, `debug`
- **fcs_screen**: `container_engine`, `workdir`, `infile`, `image`, `prok`, `fcs_script`, `outfile`, `debug`
- **fcs_gx_purge**: `workdir`, `db`, `input`, `taxid`, `outfile`, `debug`
- **sourpurge**: `workdir`, `cpus`, `read1`, `read2`, `sourdb`, `sourdb_type`, `input`, `kmer`, `phylum`, `mincovpct`, `outfile`, `taxonomy`, `debug`
- **rmdup**: `workdir`, `cpus`, `input`, `percent_id`, `percent_cov`, `minlen`, `exhaustive`, `debug`, `out`
- **polish**: `method`, `memory`, `cpus`, `read1`, `read2`, `longreads`, `workdir`, `infile`, `outfile`, `debug`. Each `--method` (`pypolca`, `polypolish`, `nextpolish2`, `racon`) is isolated into its own `run_<method>()` function in `aaftf/polish.py`.
- **sort**: `input`, `minlen`, `out`, `name`
- **assess**: `input`, `report`, `telomere_monomer`, `telomere_n_repeat`, `telomere_window`
- **depth**: `input`, `out`, `read1`, `read2`, `longreads`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `min_contig_len`, `no_plot`, `plot_format`
- **mito**: `workdir`, `read1`, `read2`, `seed`, `minlen`, `maxlen`, `out`, `subsample`, `memory`, `debug`
- **fix_tbl**: `table`, `report`, `output`, `debug`

## Pipeline step defaults

`pipeline.py` runs trim → filter → assemble → vecscreen → sourpurge → rmdup → sort → assess (mito, fcs_screen, fcs_gx_purge, polish and depth are optional and not run; polish is not recommended for short-read-only assemblies). Each step is called with `_step_kwargs(name, shared, **step_options)`:
- it starts from that step's CLI defaults, read from its `_menu.py` parser (`_subcommand_defaults()`), so every `run()` parameter is present and the pipeline never drifts from `AAFTF <step>` defaults — do not hard-code step settings in `pipeline.py`
- `shared` pipeline options (`cpus`, `memory`, `workdir`, `debug`, `quiet`) go only to steps that have that option; a `None` value (e.g. no `--memory`) keeps the step default
- `step_options` are the pipeline options specific to that step plus the file names chaining the steps
- `tests/test_pipeline.py` checks that each step receives exactly its CLI options, with defaults kept unless a pipeline option overrides them

## Pipeline Workflow

Recommended order for full assembly QC:

```
trim → [mito optional, PE only] → filter → assemble → vecscreen
  → [fcs_screen / fcs_gx_purge optional] → sourpurge → rmdup
  → [polish optional, long reads] → sort → assess → [depth optional]
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
| polish (optional) | rmdup FASTA | `{base}.polish.fasta` |
| sort | rmdup (or polished) FASTA | `{base}.final.fasta` |
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

**New file:** `aaftf/depth.py`; registered in `main.py`.

**Workflow:**
1. Counts reads in each input FASTQ
2. Maps Illumina reads with `minimap2 -ax sr` (default) or `bwa mem` (`--aligner bwa`); maps long reads with `minimap2 -ax map-ont` (default) / `map-pb` / `map-hifi` (`--longread_type`)
3. Merges BAMs when both read types are provided (`samtools merge`)
4. Runs `samtools flagstat` per read type, then `mosdepth` on the combined BAM
5. Parses `mosdepth.summary.txt` for per-contig mean depths and `mosdepth.global.dist.txt` for coverage breadth
6. Flags contigs with mean depth > assembly_mean + 3×SD as possible contaminants/organelles (uses population SD; contigs are the full assembly population, not a sample)

**Key `args` attributes:** `input`, `out`, `read1`, `read2`, `longreads`, `longread_preset`, `aligner`, `cpus`, `workdir`, `debug`, `min_contig_len`, `no_plot`, `plot_format`

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

### CLI consistency fixes (main.py)

- Added `--pipe` and `-v/--debug` to `fix_tbl` parser
- Added `--pipe` to `pipeline` parser
- Removed dead `runall` stub from `run_subtool()`
- Fixed assembler help text: "Default: unicycler" → "Default: spades"

### pipeline.py refactor

- `sort_args` and `assess_args` now use `create_namespace()` (consistent with all other steps)
- Assembly-specific attributes (`merged`, `isolate`, `careful`) now passed via `required_args` before `create_namespace()` call

### Bioinformatics enhancements

- `assess.py`: `find_telomere()` window is now configurable (`--telomere_window`, default 200 bp); added N-gap count and total N bases to stats output
- `rmdup.py`: now computes and reports both N50 and N75; filtering threshold uses N75 (as before)
- `depth.py`: added `--min_contig_len` (default 500 bp) to exclude short contigs from outlier depth analysis

## Running Tests

```bash
# Unit tests only (no external tools required)
conda run -n base python -m pytest tests/ -m unit -v

# All tests (integration tests require external tool environment)
conda run -n base python -m pytest tests/ -v

# Compile-check all modules
python -m py_compile aaftf/filter.py aaftf/polish.py aaftf/trim.py \
    aaftf/vecscreen.py aaftf/sourpurge.py aaftf/fcs_screen.py \
    aaftf/depth.py aaftf/main.py aaftf/pipeline.py \
    aaftf/assess.py aaftf/rmdup.py aaftf/_menu.py
```

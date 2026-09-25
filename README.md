# AAFTF - Automatic Assembly For The Fungi
*Authors: Jason Stajich and Jon Palmer*

![AAFTF logo](docs/AAFTF.png)

# Requirements
Python >=3.10 and the external tools below; all can be installed from conda (bioconda) or with pixi (see Install).
AAFTF needs samtools >=1.13. Run `AAFTF dependency` to see which tools are found.

## read aligners for polishing and depth of coverage calculation
- bwa - https://github.com/lh3/bwa
- minimap2 - https://github.com/lh3/minimap2
- bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml (Optional; not default)
- BBTools - https://github.com/bbushnell/BBTools

## QC and trimming
- BBTools - https://bbmap.org/ - supports read-level filtering for contamination and vector/primer
- Trimmomatic - https://github.com/usadellab/Trimmomatic (Optional; not default)
- fastp - alternative (preferred) read trimming and quality control https://github.com/OpenGene/fastp

## Assemblers
- SPAdes - https://github.com/ablab/spades
- megahit - https://github.com/voutcn/megahit
- NOVOplasty - https://github.com/ndierckx/NOVOPlasty for MT genome assembly
- unicycler - https://github.com/rrwick/Unicycler (which runs spades)

## Assembly Contamination screening support
- [sourmash](https://pubmed.ncbi.nlm.nih.gov/31508216/) (>=v4)- https://sourmash.readthedocs.io/ (install via conda/pip)
- NCBI BLAST+ - https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- [ncbi-fcs](https://pubmed.ncbi.nlm.nih.gov/38409096/) (for vector screening) - https://github.com/ncbi/fcs/
- [ncbi-fcs-gx](https://pubmed.ncbi.nlm.nih.gov/38409096/) (for contaminant filtering, alternative to sourmash, requires large memory or SSD drive) https://github.com/ncbi/fcs-gx


## Assembly polishing
- [polca](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007981) via [pypolca](https://github.com/replikation/pypolca) -
  a Python reimplementation of MaSuRCA's POLCA algorithm that works with any modern samtools
- [Polypolish](https://github.com/rrwick/Polypolish) - alignment-filtering short-read polisher
- [NextPolish2](https://github.com/Nextomics/NextPolish2) - repeat-aware polishing of HiFi assemblies using a short-read k-mer (yak) database
- [Racon](https://github.com/lbcb-sci/racon) - long-read polishing

## Depth of coverage
- [mosdepth](https://github.com/brentp/mosdepth) and samtools


# Authors
* Jason Stajich [@hyphaltip](https://github.com/hyphaltip) - http://lab.stajich.org, [@hyphaltip.bsky](https://bsky.app/profile/hyphaltip.bsky.social)
* Jon Palmer [@nextgenusfs](https://github.com/nextgenusfs) - [@jonpalmer.bsky](https://bsky.app/profile/jonpalmer.bsky.social)

# Citation
Palmer JM and Stajich JE. (2023). Automatic assembly for the fungi (AAFTF): genome assembly pipeline (v0.5.0). Zenodo. doi: 10.5281/zenodo.1620526

# Install
With conda, from a checkout of this repository:

```
# only what `AAFTF pipeline` needs with its default settings; installs AAFTF from PyPI
$ conda env create -f environment.yml && conda activate aaftf

# or every tool AAFTF can use, with AAFTF installed from this checkout (editable)
$ conda env create -f environment.dev.yml && conda activate aaftf-dev
```

Or with [pixi](https://pixi.sh), which uses the locked versions in `pixi.lock`:

```
$ pixi install -e complete        # environments: default, complete, dev
$ pixi shell -e complete
```

`default` has only the default-pipeline tools, `complete` adds every optional tool, and `dev` adds
the test and lint tools. `AAFTF dependency` reports what an environment is missing. See
`docs/installation.rst` for details, including the Docker and Singularity images.

AAFTF caches its reference databases in `~/.cache/aaftf` unless `AAFTF_DB` is set. Some are several GB,
so point `AAFTF_DB` at a folder with plenty of space (it can list several folders separated by `:`, like
`$PATH`; databases are read from the first folder that has them and downloaded into the first writable one):
```
$ export AAFTF_DB=/path/with/space/aaftf_db   # or a shared, system-wide location
$ AAFTF database                                        # list databases and where they are stored
$ AAFTF database phix univec euks proks mitodb sm_gbk   # download by short name (or file name, or 'all')
```

To run ncbi-fcs or ncbi-fcs-gx in AAFTF through singularity will need to have that installed in system or environment.
The fcs gx database will need to be downloaded and requires large memory machines.
More instructions coming for simplicity of install/testing.

## Notes
This is partially a python re-write of [JAAWS](https://github.com/nextgenusfs/jaaws) which was a unix shell based cleanup and assembly tool written by Jon.

## Steps / Procedures
Setup: `dependency` (check installed tools) and `database` (list/download reference databases).

1. trim                Trim FASTQ input reads - with bbduk (default), Trimmomatic or fastp
2. mito                (Optional) De novo assemble mitochondrial genome - with NOVOPlasty
3. filter              Filter contaminanting reads - with bbduk (default), bowtie2, bwa or minimap2
4. assemble            Assemble reads - with SPAdes (default), MEGAHIT or Unicycler
5. vecscreen           Vector and Contaminant Screening of assembled contigs - with BlastN based method to replicate NCBI screening
6a. sourpurge          Purge contigs based on sourmash results - with sourmash
6b. fcs_screen         (Optional) NCBI FCS-adaptor vector screening
6c. fcs_gx_purge       (Optional) Purge contigs based on NCBI fcs-gx tool. Note this runs MUCH faster with large memory.
7. rmdup               Remove duplicate contigs - using minimap2 to find duplicates
8. polish              (Optional) Polish contig sequences - uses Polypolish (default), pypolca, NextPolish2, or Racon; recommended only with long reads
9. sort                Sort contigs by length and rename FASTA headers
10. assess             Assess completeness of genome assembly
11. depth              (Optional) Calculate read depth of coverage across assembled contigs
12. pipeline           Run trim, filter, assemble, vecscreen, sourpurge, rmdup, sort and assess in one go

Annotation: `fix_tbl` fixes .tbl feature-table offsets after contigs were trimmed.


# Typical runs


## Trimming and Filtering

Trimming options spelled out:
```
usage: AAFTF trim [-h] -1 FASTQ [-2 FASTQ] [-o PREFIX] [-ml BP] [-aq INT]
                  [--cutfront] [--cuttail] [--cutright]
                  [--method {bbduk,trimmomatic,fastp}] [-c INT] [-m GB] [-q]
                  [-v] [--trimmomatic_adaptors TRIMMOMATIC_ADAPTORS]
                  [--trimmomatic_clip TRIMMOMATIC_CLIP]
                  [--trimmomatic_leadingwindow TRIMMOMATIC_LEADINGWINDOW]
                  [--trimmomatic_trailingwindow TRIMMOMATIC_TRAILINGWINDOW]
                  [--trimmomatic_slidingwindow TRIMMOMATIC_SLIDINGWINDOW]
                  [--trimmomatic_quality TRIMMOMATIC_QUALITY] [--dedup]
                  [--merge]

This command trims reads in FASTQ format to remove low quality reads and trim
adaptor sequences

options:
  -h, --help            show this help message and exit

required arguments:
  -1 FASTQ, --read1 FASTQ
                        Read 1 (forward) FASTQ, or single-end FASTQ

optional arguments:
  -2 FASTQ, --read2 FASTQ
                        Read 2 (reverse) FASTQ for paired-end data
  -o PREFIX, --out PREFIX
                        Output file prefix; default: read 1's file name up to
                        its first '_' (or first '.' if it has none)
  -ml BP, --minlen BP   Minimum read length to keep after trimming (default:
                        75)
  -aq INT, --avgqual INT
                        Average Quality of reads must be > than this (default:
                        10)
  --cutfront            Run fastp 5' trimming based on quality. WARNING: this
                        operation will interfere deduplication for SE data
  --cuttail             Run fastp 3' trimming based on quality. WARNING: this
                        operation will interfere deduplication for SE data
  --cutright            Run fastp move a sliding window from front to tail, if
                        meet one window with mean quality < threshold.
                        WARNING: this operation will interfere deduplication
                        for SE data
  --method {bbduk,trimmomatic,fastp}
                        Trimming method (default: bbduk)
  -c INT, --cpus INT    Number of CPUs/threads to use (default: 1)
  -m GB, --memory GB    Max memory in GB (default: 8)
  -q, --quiet           Only show warnings and errors
  -v, --verbose         Show debug messages and tool stderr, and keep
                        temporary working directories

Trimmomatic options:
  --trimmomatic_adaptors TRIMMOMATIC_ADAPTORS
                        Trimmomatic adaptor file (default: TruSeq3-PE.fa)
  --trimmomatic_clip TRIMMOMATIC_CLIP
                        Trimmomatic ILLUMINACLIP argument (default: 2:30:10)
  --trimmomatic_leadingwindow TRIMMOMATIC_LEADINGWINDOW
                        Trimmomatic window processing arguments (default: 3)
  --trimmomatic_trailingwindow TRIMMOMATIC_TRAILINGWINDOW
                        Trimmomatic window processing arguments (default: 3)
  --trimmomatic_slidingwindow TRIMMOMATIC_SLIDINGWINDOW
                        Trimmomatic window processing arguments (default:
                        4:15)
  --trimmomatic_quality TRIMMOMATIC_QUALITY
                        Trimmomatic quality encoding -phred33 or phred64
                        (default: phred33)

Fastp options:
  --dedup               Run fastp deuplication of fastq reads (default uses
                        ~4gb mem)
  --merge               Merge paired end reads
```

Example usage:
```
MEM=128 # 128gb
BASE=STRAINX
READSDIR=reads
TRIMREAD=reads_trimmed
CPU=8
AAFTF trim --method bbduk --memory $MEM -c $CPU \
 --read1 $READSDIR/${BASE}_R1.fq.gz --read2 $READSDIR/${BASE}_R2.fq.gz \
  -o $TRIMREAD/${BASE}
# this step make take a lot of memory depending on how many filtering libraries you use
AAFTF filter -c $CPU --memory $MEM --aligner bbduk \
	  -o $TRIMREAD/${BASE} --read1 $TRIMREAD/${BASE}_1P.fastq.gz --read2 $TRIMREAD/${BASE}_2P.fastq.gz
```

## Assembly

The specified assembler can be made through the `--method` option.
The full set of options are below.

```
usage: AAFTF assemble [-h] -1 FASTQ -o FASTA [-2 FASTQ] [-w DIR]
                      [--method {spades,megahit,unicycler}] [--merged MERGED]
                      [--tmpdir DIR] [--assembler_args ARG] [-c INT] [-m GB]
                      [-q] [-v] [--no-careful] [--no-isolate] [-lr FASTQ]

Run assembler on cleaned reads

options:
  -h, --help            show this help message and exit

required arguments:
  -1 FASTQ, --read1 FASTQ
                        Read 1 (forward) FASTQ, or single-end FASTQ
  -o FASTA, --out FASTA
                        Output assembly FASTA

optional arguments:
  -2 FASTQ, --read2 FASTQ
                        Read 2 (reverse) FASTQ for paired-end data
  -w DIR, --workdir DIR
                        Working directory for intermediate files; a temporary
                        one is created and removed afterwards (kept with -v)
                        when not given
  --method {spades,megahit,unicycler}
                        Assembly method (default: spades)
  --merged MERGED       Merged reads from flash or fastp or just single end
                        reads
  --tmpdir DIR          Temporary directory for the assembler
  --assembler_args ARG  Extra argument passed to the assembler (repeat for
                        several)
  -c INT, --cpus INT    Number of CPUs/threads to use (default: 1)
  -m GB, --memory GB    Max memory in GB for the assembler (default: 32)
  -q, --quiet           Only show warnings and errors
  -v, --verbose         Show debug messages and tool stderr, and keep
                        temporary working directories

SPAdes options:
  --no-careful          Disable --careful mode in spades (Default: --careful
                        is on) (default: True)
  --no-isolate          Disable --isolate mode in spades (Default: --isolate
                        is on) (default: True)

Unicycler options:
  -lr FASTQ, --longreads FASTQ
                        Long-read FASTQ (PacBio or ONT)
```

```
CPU=24
MEM=96
READ1=$TRIMREAD/${BASE}_filtered_1.fastq.gz
READ2=$TRIMREAD/${BASE}_filtered_2.fastq.gz
WORKDIR=working_AAFTF
OUTDIR=genomes
ASMFILE=$OUTDIR/${BASE}.spades.fasta
mkdir -p $WORKDIR $OUTDIR
AAFTF assemble -c $CPU --memory $MEM \
	  --read1 $READ1 --read2 $READ2  \
	   -o $ASMFILE -w $WORKDIR/spades_$BASE
```

## vectrim

```
CPU=16
MEM=16
READ1=$TRIMREAD/${BASE}_filtered_1.fastq.gz
READ2=$TRIMREAD/${BASE}_filtered_2.fastq.gz
WORKDIR=working_AAFTF
OUTDIR=genomes
ASMFILE=$OUTDIR/${BASE}.spades.fasta
VECTRIM=$OUTDIR/${BASE}.vecscreen.fasta
mkdir -p $WORKDIR $OUTDIR
AAFTF vecscreen -c $CPU -i $ASMFILE -o $VECTRIM
```

## Depth of Coverage

The `depth` subtool maps reads to the final assembly and
computes per-contig depth statistics using mosdepth.  It requires `samtools` and `mosdepth`
plus at least one of `minimap2` (default for both Illumina and long reads) or `bwa`.

```
AAFTF depth -i genome.final.fasta \
    --read1 reads_1P.fastq.gz --read2 reads_2P.fastq.gz \
    -c $CPU -o coverage_report.txt
```

Long reads can be added alongside or instead of Illumina reads:

```
AAFTF depth -i genome.final.fasta \
    --read1 reads_1P.fastq.gz --read2 reads_2P.fastq.gz \
    --longreads nanopore.fastq.gz \
    -c $CPU -o coverage_report.txt
```

The report (`coverage_stats.txt` by default) contains three sections:

1. **Read input summary** — per-file read counts and `samtools flagstat` alignment rates
2. **Whole-assembly coverage** — two mean depth estimates:
   - *mosdepth global* (length-weighted: total bases covered ÷ assembly length)
   - *per-contig arithmetic mean* (unweighted mean of per-contig means)
   — plus percentage of bases covered at ≥1x
3. **Per-contig depth table** (sorted by depth, descending) — each contig is flagged:
   - `** OUTLIER` — mean depth > assembly mean + 3 SD (likely contaminant or organelle)
   - `ELEVATED`  — mean depth between 2 SD and 3 SD above mean (worth inspecting)

When matplotlib is available, three coverage plots are also produced alongside the report.

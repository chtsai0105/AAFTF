"""Consolidated argparse subcommand parser (menu) definitions for AAFTF.

Each ``<name>_menu(subparsers)`` function below builds and registers one
AAFTF subcommand's argparse parser via ``subparsers.add_parser(...)``, and
binds that subtool's ``run()`` function to it via ``parser_x.set_defaults(
func=<module>.run)``. This keeps all CLI surface/wiring in one place,
separate from each subcommand module's own ``run(**kwargs)`` execution logic.
``AAFTF_main.py`` invokes the selected subtool via ``args.func(**vars(args))``.

``menu_common_args()`` adds the three arguments common to every subcommand
(``-v/--debug``, ``--pipe``, ``-q/--quiet``). Call it last, passing the
subcommand's own "optional arguments" group (not the parser itself), so these
common flags render as the final entries in that group instead of appearing
in a separate leading section.
"""

import argparse as ap

import AAFTF.assemble as assemble
import AAFTF.assess as assess
import AAFTF.check_dependencies as check_dependencies
import AAFTF.depth as depth
import AAFTF.download as download
import AAFTF.fcs_gx_purge as fcs_gx_purge
import AAFTF.fcs_screen as fcs_screen
import AAFTF.filter as aaftf_filter
import AAFTF.fix_tbl as fix_tbl
import AAFTF.mito as mito
import AAFTF.pipeline as pipeline
import AAFTF.polish as polish
import AAFTF.rmdup as rmdup
import AAFTF.sort as aaftf_sort
import AAFTF.sourpurge as sourpurge
import AAFTF.trim as trim
import AAFTF.vecscreen as vecscreen
from AAFTF.utility import CustomHelpFormatter


def menu_common_args(target):
    """Add the arguments common to every AAFTF subcommand parser.

    ``target`` is normally a subcommand's "optional arguments" group (from
    ``parser.add_argument_group("optional arguments")``). Call this only
    after all of that subcommand's own optional arguments have been added,
    so ``-v/--debug``, ``--pipe``, and ``-q/--quiet`` render as the final
    entries of that group.
    """
    target.add_argument("--pipe", action="store_true", help="AAFTF is running in pipeline mode")
    target.add_argument("-q", "--quiet", action="store_true", dest="quiet", help="Do not output warnings to stderr")
    target.add_argument("-v", "--debug", action="store_true", help="Provide debugging messages")
    return target


def download_menu(subparsers):
    """Add the download subcommand parser."""
    parser_download = subparsers.add_parser(
        "download",
        formatter_class=CustomHelpFormatter,
        description="Download reference databases to a persistent directory.",
        help="Download AAFTF reference databases",
    )
    optional = parser_download.add_argument_group("optional arguments")
    optional.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist",
    )
    optional.add_argument(
        "--skip-core",
        action="store_true",
        help="Skip downloading core contamination databases (UniVec, PhiX, contaminants, mitochondria)",
    )
    optional.add_argument(
        "--skip-sourmash",
        action="store_true",
        help="Skip downloading sourmash taxonomy databases",
    )
    optional.add_argument(
        "--sourdb-type",
        type=str,
        default="gbk",
        choices=["gbk", "gtdb", "gtdbrep", "all"],
        dest="sourdb_type",
        help="Which sourmash database(s) to download",
    )
    optional.add_argument(
        "--skip-fcs",
        action="store_true",
        help="Skip downloading NCBI FCS-adaptor resources",
    )
    optional.add_argument(
        "--list",
        action="store_true",
        dest="list_db",
        help="List database files already present in AAFTF_DB (with sizes) instead of downloading",
    )
    menu_common_args(optional)
    parser_download.set_defaults(func=download.run)
    return parser_download


def trim_menu(subparsers):
    """Add the trim subcommand parser."""
    parser_trim = subparsers.add_parser(
        "trim",
        formatter_class=CustomHelpFormatter,
        description="This command trims reads in FASTQ format to remove low quality reads and trim adaptor sequences",
        help="Trim FASTQ input reads",
    )

    required = parser_trim.add_argument_group("required arguments")
    optional = parser_trim.add_argument_group("optional arguments")

    required.add_argument(
        "-l",
        "--left",
        type=str,
        metavar="FASTQ",
        required=True,
        help="left/forward reads of paired-end FASTQ or single-end FASTQ.",
    )

    optional.add_argument("-r", "--right", type=str, metavar="FASTQ", help="right/reverse reads of paired-end FASTQ.")

    optional.add_argument(
        "-o",
        "--out",
        type=str,
        dest="basename",
        help="Output basename, default to base name of --left reads",
    )

    optional.add_argument("-ml", "--minlen", type=int, metavar="INT", default=75, help="Minimum read length after trimming")

    optional.add_argument(
        "-aq",
        "--avgqual",
        type=int,
        metavar="INT",
        default=10,
        help="Average Quality of reads must be > than this",
    )

    optional.add_argument(
        "--cutfront",
        action="store_true",
        help="Run fastp 5' trimming based on quality. WARNING: this operation will interfere deduplication for SE data",
    )

    optional.add_argument(
        "--cuttail",
        action="store_true",
        help="Run fastp 3' trimming based on quality. WARNING: this operation will interfere deduplication for SE data",
    )

    optional.add_argument(
        "--cutright",
        action="store_true",
        help="Run fastp move a sliding window from front to tail, if meet one window with mean quality < threshold. \nWARNING: this operation will interfere deduplication for SE data",
    )

    optional.add_argument("--method", default="bbduk", choices=["bbduk", "trimmomatic", "fastp"], help="Program to use for adapter trimming")

    optional.add_argument("-c", "--cpus", type=int, metavar="cpus", default=1, help="Number of CPUs/threads to use.")
    optional.add_argument("-m", "--memory", type=int, dest="memory", help="Max Memory (in GB)")
    menu_common_args(optional)

    trimmomatic_group = parser_trim.add_argument_group(title="Trimmomatic options")

    trimmomatic_group.add_argument("--trimmomatic_adaptors", default="TruSeq3-PE.fa", help="Trimmomatic adaptor file")

    trimmomatic_group.add_argument("--trimmomatic_clip", type=str, default="2:30:10", help="Trimmomatic ILLUMINACLIP argument")

    trimmomatic_group.add_argument("--trimmomatic_leadingwindow", type=int, default="3", help="Trimmomatic window processing arguments")

    trimmomatic_group.add_argument("--trimmomatic_trailingwindow", type=int, default="3", help="Trimmomatic window processing arguments")

    trimmomatic_group.add_argument(
        "--trimmomatic_slidingwindow",
        type=str,
        default="4:15",
        help="Trimmomatic window processing arguments",
    )

    trimmomatic_group.add_argument("--trimmomatic_quality", default="phred33", help="Trimmomatic quality encoding -phred33 or phred64")

    fastp_group = parser_trim.add_argument_group(title="Fastp options")

    fastp_group.add_argument(
        "--dedup",
        action="store_true",
        help="Run fastp deuplication of fastq reads (default uses ~4gb mem)",
    )

    fastp_group.add_argument("--merge", action="store_true", help="Merge paired end reads")

    parser_trim.set_defaults(func=trim.run)
    return parser_trim


def mito_menu(subparsers):
    """Add the mito subcommand parser."""
    parser_mito = subparsers.add_parser(
        "mito",
        description="De novo assembly of mitochondrial genome using NOVOplasty, takes PE Illumina adapter trimmed data.",
        help="De novo assembly of mitochondrial genome",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_mito.add_argument_group("required arguments")
    optional = parser_mito.add_argument_group("optional arguments")

    required.add_argument("-l", "--left", required=True, help="Left (Forward) reads")

    required.add_argument("-r", "--right", required=True, help="Right (Reverse) reads")

    optional.add_argument("-o", "--out", type=str, default="mito.fasta", help="Output FASTA file for mitochondrial genome")

    optional.add_argument(
        "-w",
        "--workdir",
        "--tmpdir",
        type=str,
        dest="workdir",
        help="Temporary directory to store datafiles and processes in",
    )

    optional.add_argument("--minlen", default=10000, type=int, help="Minimum expected genome size")

    optional.add_argument("--maxlen", default=100000, type=int, help="Maximum expected genome size")

    optional.add_argument("-s", "--seed", help="Seed sequence, ie related mitochondrial genome. default: A. nidulans")

    optional.add_argument("--starting", help="FASTA file of start sequence, rotate genome to, default COB")

    optional.add_argument("--reference", help="Run NOVOplasty in reference mode")

    menu_common_args(optional)

    parser_mito.set_defaults(func=mito.run)
    return parser_mito


def filter_menu(subparsers):
    """Add the filter subcommand parser."""
    parser_filter = subparsers.add_parser(
        "filter",
        description="Filter reads which match contaminant databases such as phiX",
        help="Filter contaminanting reads",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_filter.add_argument_group("required arguments")
    optional = parser_filter.add_argument_group("optional arguments")

    required.add_argument("-l", "--left", required=True, help="Left (Forward) reads")

    optional.add_argument("-r", "--right", help="Right (Reverse) reads")

    optional.add_argument("-o", "--out", dest="basename", type=str, help="Output basename")

    optional.add_argument(
        "-w",
        "--workdir",
        "--tmpdir",
        type=str,
        dest="workdir",
        help="Temporary directory to store datafiles and processes in",
    )

    optional.add_argument(
        "--aligner",
        default="bbduk",
        choices=["bbduk", "bowtie2", "bwa", "minimap2"],
        help="Aligner to use to map reads to contamination database",
    )

    optional.add_argument("-a", "--screen_accessions", type=str, nargs="*", help="Genbank accession number(s) to screen out from initial reads.")

    optional.add_argument("-u", "--screen_urls", type=str, nargs="*", help="URLs to download and screen out initial reads.")

    optional.add_argument("-s", "--screen_local", type=str, nargs="+", help="Local FASTA file(s) to use contamination screen")

    optional.add_argument("-c", "--cpus", type=int, metavar="cpus", default=1, help="Number of CPUs/threads to use.")

    optional.add_argument("-m", "--memory", type=int, dest="memory", help="Max Memory (in GB)")

    menu_common_args(optional)

    parser_filter.set_defaults(func=aaftf_filter.run)
    return parser_filter


def assemble_menu(subparsers):
    """Add the assemble subcommand parser."""
    parser_asm = subparsers.add_parser(
        "assemble",
        description="Run assembler on cleaned reads",
        help="Assemble reads",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_asm.add_argument_group("required arguments")
    optional = parser_asm.add_argument_group("optional arguments")

    required.add_argument(
        "-l",
        "--left",
        required=True,  # every implemented --method (spades/dipspades/megahit/unicycler) requires this
        help="Left (Forward) reads",
    )

    required.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,  # think about sensible replacement in future
        help="Output assembly FASTA",
    )

    optional.add_argument("-r", "--right", help="Right (Reverse) reads")

    optional.add_argument("-w", "--workdir", type=str, dest="workdir", help="assembly output directory")

    optional.add_argument(
        "--method",
        type=str,
        choices=["spades", "dipspades", "megahit", "unicycler"],
        default="spades",
        help="Assembly method: spades, dipspades, megahit, unicycler",
    )

    optional.add_argument("--merged", dest="merged", help="Merged reads from flash or fastp or just single end reads")
    optional.add_argument("--tmpdir", type=str, help="Assembler temporary dir")
    optional.add_argument("--assembler_args", action="append", help="Additional SPAdes/Megahit arguments")
    optional.add_argument("-c", "--cpus", type=int, metavar="cpus", default=1, help="Number of CPUs/threads to use.")
    optional.add_argument(
        "-m",
        "--memory",
        type=str,
        dest="memory",
        default="32",
        help="Memory (in GB) setting for SPAdes",
    )

    menu_common_args(optional)

    spades_group = parser_asm.add_argument_group(title="SPAdes options")

    spades_group.add_argument(
        "--no-careful",
        action="store_false",
        default=True,
        dest="careful",
        help="Disable --careful mode in spades (Default: --careful is on)",
    )

    spades_group.add_argument(
        "--no-isolate",
        action="store_false",
        default=True,
        dest="isolate",
        help="Disable --isolate mode in spades (Default: --isolate is on)",
    )

    dipspades_group = parser_asm.add_argument_group(title="dipSPAdes options")

    dipspades_group.add_argument("--haplocontigs", dest="haplocontigs", default=False, action="store_true", help="For dipSPAdes take the haplocontigs file")

    unicycler_group = parser_asm.add_argument_group(title="Unicycler options")

    unicycler_group.add_argument("-lr", "--longreads", help="Long Read fastq (pacbio or ONT)")

    parser_asm.set_defaults(func=assemble.run)
    return parser_asm


def vecscreen_menu(subparsers):
    """Add the vecscreen subcommand parser."""
    parser_vecscreen = subparsers.add_parser(
        "vecscreen",
        description="Screen contigs for vector and common contaminantion",
        help="BLASTN Vector and Contaminant Screening of contigs",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_vecscreen.add_argument_group("required arguments")
    optional = parser_vecscreen.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", "--infile", type=str, required=True, dest="infile", help="Input contigs or scaffold assembly")

    required.add_argument("-o", "--outfile", type=str, required=True, help="Output vector screened and cleaned assembly")

    optional.add_argument(
        "-w",
        "--workdir",
        "--tmpdir",
        type=str,
        dest="workdir",
        help="Temporary directory to store datafiles and processes in",
    )

    optional.add_argument("-pid", "--percent_id", type=int, help="Percent Identity cutoff for vecscreen adaptor matches")

    optional.add_argument("-s", "--stringency", default="high", choices=["high", "low"], help="Stringency to filter VecScreen hits")

    optional.add_argument("-c", "--cpus", type=int, metavar="cpus", default=1, help="Number of CPUs/threads to use.")

    menu_common_args(optional)

    parser_vecscreen.set_defaults(func=vecscreen.run)
    return parser_vecscreen


def fcs_screen_menu(subparsers):
    """Add the fcs_screen subcommand parser."""
    parser_fcs_screen = subparsers.add_parser(
        "fcs_screen",
        description="Screen with NCBI fcs tool contigs for vector and common contaminantion",
        help="NCBI Foreign Contaminant Screening for Vector sequences in contigs",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_fcs_screen.add_argument_group("required arguments")
    optional = parser_fcs_screen.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", "--infile", type=str, required=True, dest="infile", help="Input contigs or scaffold assembly")

    required.add_argument("-o", "--outfile", type=str, required=True, help="Output vector screened and cleaned assembly")

    optional.add_argument(
        "-w",
        "--workdir",
        "--tmpdir",
        type=str,
        dest="workdir",
        help="Temporary directory to store datafiles and processes in",
    )

    optional.add_argument("--image", type=str, help="Container file (or will download and look in AAFTF_DB)")

    optional.add_argument(
        "--container_engine",
        type=str,
        default="singularity",
        choices=["singularity", "docker"],
        help="Container engine used to run fcs-adaptor",
    )

    optional.add_argument("--prok", action="store_true", help="Run in Prokaryote matching mode")

    optional.add_argument("--fcs_script", type=str, help="location of the run_fcsadaptor.sh script (or will download automatically)")

    menu_common_args(optional)

    parser_fcs_screen.set_defaults(func=fcs_screen.run)
    return parser_fcs_screen


def fcs_gx_purge_menu(subparsers):
    """Add the fcs_gx_purge subcommand parser."""
    parser_fcsgx = subparsers.add_parser(
        "fcs_gx_purge",
        description="Purge contigs based on fcs_gx results",
        help="Purge contigs based on contamination search with fcs_gx",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_fcsgx.add_argument_group("required arguments")
    optional = parser_fcsgx.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", type=str, required=True, help="Input contigs or scaffold assembly")

    required.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,  # think about sensible replacement in future
        help="Output fcs_gx cleaned assembly",
    )

    optional.add_argument(
        "-w",
        "--workdir",
        "--tmpdir",
        type=str,
        dest="workdir",
        help="Temporary directory to store datafiles and processes in",
    )

    optional.add_argument(
        "-t",
        "--taxid",
        type=int,
        default=4890,
        help="NCBI Taxonomy ID for contamaination matches, i.e. 4890 for Ascomycota",
    )

    optional.add_argument("-d", "--db", default="/my_tmpfs/gxdb/all", help="gxdb database path")

    menu_common_args(optional)

    parser_fcsgx.set_defaults(func=fcs_gx_purge.run)
    return parser_fcsgx


def sourpurge_menu(subparsers):
    """Add the sourpurge subcommand parser."""
    parser_sour = subparsers.add_parser(
        "sourpurge",
        description="Purge contigs based on sourmash results",
        help="Purge contigs based on sourmash results",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_sour.add_argument_group("required arguments")
    optional = parser_sour.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", type=str, required=True, help="Input contigs or scaffold assembly")

    required.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,  # think about sensible replacement in future
        help="Output sourmash cleaned assembly",
    )

    required.add_argument("-p", "--phylum", required=True, nargs="+", help="Phylum or Phyla to keep matches, i.e. Ascomycota")

    optional.add_argument(
        "-w",
        "--workdir",
        "--tmpdir",
        type=str,
        dest="workdir",
        help="Temporary directory to store datafiles and processes in",
    )

    optional.add_argument("-l", "--left", help="Left (Forward) reads")

    optional.add_argument("-r", "--right", help="Right (Reverse) reads")

    optional.add_argument("--sourdb", help="SourMash LCA taxonomy database (defaults to k-31)")

    optional.add_argument("-k", "--kmer", default="31", help="SourMash LCA kmersize when taxonomy database was built")

    optional.add_argument("-mc", "--mincovpct", default=5, type=int, help="Minimum percent of N50 coverage to remove")

    optional.add_argument(
        "--sourdb_type",
        default="gbk",
        choices=["gbk", "gtdbrep", "gtdb"],
        help="Which sourpurge database to use.",
    )

    optional.add_argument("--just-show-taxonomy", dest="taxonomy", action="store_true", help="Show taxonomy information and exit")
    optional.add_argument("-c", "--cpus", type=int, metavar="cpus", default=1, help="Number of CPUs/threads to use.")

    menu_common_args(optional)

    parser_sour.set_defaults(func=sourpurge.run)
    return parser_sour


def rmdup_menu(subparsers):
    """Add the rmdup subcommand parser."""
    parser_rmdup = subparsers.add_parser(
        "rmdup",
        description="Remove duplicate contigs",
        help="Remove duplicate contigs",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_rmdup.add_argument_group("required arguments")
    optional = parser_rmdup.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", type=str, required=True, help="Input Assembly fasta file(contigs or scaffolds)")

    required.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="Output new version of assembly with duplicated contigs/scaffolds removed",
    )

    optional.add_argument("-c", "--cpus", type=int, metavar="cpus", default=1, help="Number of CPUs/threads to use.")

    optional.add_argument(
        "-w",
        "--workdir",
        "--tmpdir",
        type=str,
        dest="workdir",
        help="Temporary directory to store datafiles and processes in",
    )

    optional.add_argument(
        "-pid",
        "--percent_id",
        type=int,
        dest="percent_id",
        default=95,
        help="Percent Identity used in matching contigs for redundancy",
    )

    optional.add_argument(
        "-pcov",
        "--percent_cov",
        type=int,
        dest="percent_cov",
        default=95,
        help="Coverage of contig used to decide if it is redundant",
    )

    optional.add_argument("-ml", "--minlen", type=int, default=500, help="Minimum contig length to keep, shorter ones are dropped")

    optional.add_argument(
        "--exhaustive",
        action="store_true",
        help="Compute overlaps for every contig, otherwise only process contigs for L75 and below",
    )

    menu_common_args(optional)

    parser_rmdup.set_defaults(func=rmdup.run)
    return parser_rmdup


def polish_menu(subparsers):
    """Add the polish subcommand parser."""
    parser_polish = subparsers.add_parser(
        "polish",
        description="Polish contig sequences with pypolca, Polypolish, or NextPolish2",
        help="Polish contig sequences with short and/or long reads",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_polish.add_argument_group("required arguments")

    required.add_argument("-i", "--infile", "--input", type=str, dest="infile", required=True, help="Input contigs or scaffold assembly")

    shortread_group = parser_polish.add_argument_group(title="polypolish / pypolca / NextPolish2 required arguments")

    shortread_group.add_argument("-l", "--left", help="Left (Forward) reads; required for short read polishing methods")

    shortread_group.add_argument("-r", "--right", help="Right (Reverse) reads; required for short read polishing methods")

    longread_group = parser_polish.add_argument_group(title="NextPolish2 / Racon required arguments")

    longread_group.add_argument("-lr", "--longreads", help="Long Read FASTQ (PacBio or ONT/HiFi); required for NextPolish2 and Racon")

    optional = parser_polish.add_argument_group("optional arguments")

    optional.add_argument("-o", "--out", "--outfile", type=str, dest="outfile", help="Output a Polished assembly")

    optional.add_argument(
        "-w",
        "--workdir",
        "--tmpdir",
        type=str,
        dest="workdir",
        help="Temporary directory to store datafiles and processes in",
    )

    optional.add_argument(
        "--method",
        type=str,
        choices=["polypolish", "pypolca", "nextpolish2", "racon"],
        default="polypolish",
        help="Polishing method: polypolish, pypolca, nextpolish2, racon",
    )

    optional.add_argument("-c", "--cpus", type=int, metavar="cpus", default=1, help="Number of CPUs/threads to use.")

    menu_common_args(optional)

    pypolca_group = parser_polish.add_argument_group(title="pypolca options")

    pypolca_group.add_argument("-m", "--memory", type=int, default=16, dest="memory", help="Max Memory (in GB)")

    parser_polish.set_defaults(func=polish.run)
    return parser_polish


def sort_menu(subparsers):
    """Add the sort subcommand parser."""
    parser_sort = subparsers.add_parser(
        "sort",
        description="Sort contigs by length and rename FASTA headers",
        help="Sort contigs by length and rename FASTA headers",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_sort.add_argument_group("required arguments")
    optional = parser_sort.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", "--infile", required=True, dest="input", help="Input genome assembly FASTA")

    required.add_argument("-o", "--out", "--output", required=True, dest="out", help="Output genome assembly FASTA")

    optional.add_argument("-ml", "--minlen", type=int, default=0, help="Minimum contig length to keep, shorter ones are dropped")

    optional.add_argument("-n", "--name", "--basename", default="scaffold", dest="name", help="Basename to rename FASTA headers")

    menu_common_args(optional)

    parser_sort.set_defaults(func=aaftf_sort.run)
    return parser_sort


def assess_menu(subparsers):
    """Add the assess subcommand parser."""
    parser_assess = subparsers.add_parser(
        "assess",
        description="Assess completeness of genome assembly",
        help="Assess completeness of genome assembly",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_assess.add_argument_group("required arguments")
    optional = parser_assess.add_argument_group("optional arguments")

    required.add_argument(
        "-i",
        "--input",
        "--infile",
        required=True,
        help="Input genome assembly to test completeness and provide summary statistics",
    )

    optional.add_argument("-r", "--report", type=str, help="Filename to save report information otherwise will print to stdout")

    optional.add_argument("-t", "--telomere_monomer", type=str, help="Telomere monomer pattern to search for.", default="TAAC{3,5}")

    optional.add_argument("-n", "--telomere_n_repeat", type=int, default=2, help="Telomere minimum number of monomer repeats.")

    optional.add_argument("--telomere_window", type=int, default=200, help="Number of bp to scan at each end for telomere repeats.")

    menu_common_args(optional)

    parser_assess.set_defaults(func=assess.run)
    return parser_assess


def fix_tbl_menu(subparsers):
    """Add the fix_tbl subcommand parser."""
    parser_fix = subparsers.add_parser(
        "fix_tbl",
        description="Fix NCBI tbl file from a trim report from NCBI-FCS",
        help="Fix the TBL file offsets from trimmimg",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_fix.add_argument_group("required arguments")
    optional = parser_fix.add_argument_group("optional arguments")

    required.add_argument("-t", "--table", "--infile", type=ap.FileType("rt"), required=True, help="Table format of annotation (NCBL tbl)")

    required.add_argument("-r", "--report", type=ap.FileType("rt"), required=True, help="FCS report (5 column CSV)")

    required.add_argument("-o", "--output", type=ap.FileType("wt"), required=True, help="Write fixed TBL file")

    menu_common_args(optional)

    parser_fix.set_defaults(func=fix_tbl.run)
    return parser_fix


def depth_menu(subparsers):
    """Add the depth subcommand parser."""
    parser_depth = subparsers.add_parser(
        "depth",
        description=("Calculate depth of coverage by mapping Illumina and/or long reads to a genome assembly with minimap2 (or bwa), then running mosdepth to compute per-contig depth statistics.  Contigs with mean depth > assembly_mean + 3*SD are flagged as possible contaminants or organellar sequences."),
        help="Calculate read depth of coverage for genome assembly",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_depth.add_argument_group("required arguments")
    optional = parser_depth.add_argument_group("optional arguments")

    required.add_argument(
        "-i",
        "--input",
        "--infile",
        required=True,
        dest="input",
        help="Input genome assembly FASTA (e.g. *.sorted.fasta)",
    )

    optional.add_argument("-c", "--cpus", type=int, metavar="cpus", default=1, help="Number of CPUs/threads to use.")

    optional.add_argument(
        "-w",
        "--workdir",
        "--tmpdir",
        type=str,
        dest="workdir",
        help="Temporary directory to store datafiles and processes in",
    )

    optional.add_argument(
        "-o",
        "--out",
        "--report",
        type=str,
        default="coverage_stats.txt",
        dest="out",
        help="Output coverage report file",
    )

    optional.add_argument(
        "-l",
        "--left",
        help="Left (Forward) Illumina reads FASTQ",
    )

    optional.add_argument(
        "-r",
        "--right",
        help="Right (Reverse) Illumina reads FASTQ",
    )

    optional.add_argument(
        "-lr",
        "--longreads",
        help="Long reads FASTQ (PacBio or ONT)",
    )

    optional.add_argument(
        "--longread_type",
        default="map-ont",
        choices=["map-ont", "map-pb", "map-hifi"],
        dest="longread_preset",
        help="minimap2 preset for long reads",
    )

    optional.add_argument(
        "--illumina_preset",
        default="sr",
        choices=["sr", "short"],
        dest="illumina_preset",
        help="minimap2 preset for Illumina reads",
    )

    optional.add_argument(
        "--aligner",
        default="minimap2",
        choices=["minimap2", "bwa"],
        help="Aligner to use for Illumina reads",
    )

    optional.add_argument(
        "--min_contig_len",
        type=int,
        default=500,
        help="Minimum contig length to include in depth outlier analysis",
    )

    optional.add_argument(
        "--plot-format",
        choices=["pdf", "svg", "png"],
        default="pdf",
        dest="plot_format",
        help="Output format for coverage plots (requires matplotlib)",
    )

    optional.add_argument(
        "--no-plot",
        action="store_true",
        dest="no_plot",
        help="Disable coverage plot generation",
    )

    optional.add_argument(
        "--quantize",
        default="0:1:4:100:200:",
        dest="quantize",
        help="mosdepth quantize bin boundaries, colon-separated with trailing colon",
    )

    optional.add_argument(
        "--quantize-labels",
        default=None,
        dest="quantize_labels",
        metavar="LABELS",
        help="Comma-separated labels for quantize bins (default: NO_COVERAGE,LOW_COVERAGE,CALLABLE,HIGH_COVERAGE,VERY_HIGH_COVERAGE for the default bins)",
    )

    menu_common_args(optional)

    parser_depth.set_defaults(func=depth.run)
    return parser_depth


def pipeline_menu(subparsers):
    """Add the pipeline subcommand parser."""
    parser_pipeline = subparsers.add_parser(
        "pipeline",
        description="Run entire AAFTF pipeline automagically",
        help="Run AAFTF pipeline",
        formatter_class=CustomHelpFormatter,
    )

    required = parser_pipeline.add_argument_group("required arguments")
    optional = parser_pipeline.add_argument_group("optional arguments")

    required.add_argument("-l", "--left", type=str, required=True, help="left/forward reads of paired-end FASTQ or single-end FASTQ.")

    required.add_argument("-o", "--out", type=str, required=True, dest="basename", help="Output basename, default to base name of --left reads")

    required.add_argument("-p", "--phylum", required=True, nargs="+", help="Phylum or Phyla to keep matches, i.e. Ascomycota")

    optional.add_argument("-c", "--cpus", type=int, metavar="cpus", default=1, help="Number of CPUs/threads to use.")

    optional.add_argument("--tmpdir", type=str, help="Assembler temporary dir")
    optional.add_argument("--assembler_args", action="append", help="Additional SPAdes/Megahit arguments")
    optional.add_argument("--method", type=str, default="spades", help="Assembly method: spades, dipspades, megahit")

    optional.add_argument("-r", "--right", type=str, help="right/reverse reads of paired-end FASTQ.")

    optional.add_argument("-m", "--memory", type=str, dest="memory", help="Memory (in GB) setting for SPAdes. Default is Auto")

    optional.add_argument("-ml", "--minlen", type=int, default=75, help="Minimum read length after trimming")

    optional.add_argument("-a", "--screen_accessions", type=str, nargs="*", help="Genbank accession number(s) to screen out from initial reads.")

    optional.add_argument("-u", "--screen_urls", type=str, nargs="*", help="URLs to download and screen out initial reads.")

    optional.add_argument("-mc", "--mincontiglen", type=int, default=500, help="Minimum length of contigs to keep")

    optional.add_argument("-w", "--workdir", type=str, help="temp directory")

    optional.add_argument("--sourdb", help="SourMash LCA k-31 taxonomy database")

    optional.add_argument("--mincovpct", default=5, type=int, help="Minimum percent of N50 coverage to remove")

    menu_common_args(optional)

    parser_pipeline.set_defaults(func=pipeline.run)
    return parser_pipeline


def check_dependencies_menu(subparsers):
    """Add the check_dependencies subcommand parser."""
    parser_check_deps = subparsers.add_parser(
        "check_dependencies",
        formatter_class=CustomHelpFormatter,
        description="Check whether all external tool and Python package dependencies required by AAFTF are installed.",
        help="Check that AAFTF dependencies are installed",
    )
    parser_check_deps.set_defaults(func=check_dependencies.run)
    return parser_check_deps


SUBCOMMAND_REGISTRARS = [
    download_menu,
    trim_menu,
    mito_menu,
    filter_menu,
    assemble_menu,
    vecscreen_menu,
    fcs_screen_menu,
    fcs_gx_purge_menu,
    sourpurge_menu,
    rmdup_menu,
    polish_menu,
    sort_menu,
    assess_menu,
    fix_tbl_menu,
    depth_menu,
    pipeline_menu,
    check_dependencies_menu,
]

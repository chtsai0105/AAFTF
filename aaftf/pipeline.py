"""Support pipelining of AAFTF to simplify all-in-on runs."""

import logging
import sys
from argparse import Namespace

import aaftf.assemble as assemble
import aaftf.assess as assess
import aaftf.filter as aaftf_filter
import aaftf.mito as mito
import aaftf.polish as polish
import aaftf.rmdup as rmdup
import aaftf.sort as aaftf_sort
import aaftf.sourpurge as sourpurge
import aaftf.trim as trim
import aaftf.vecscreen as vecscreen
from aaftf.utility import check_file, get_ram

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(
    left,
    basename,
    phylum,
    right=None,
    cpus=1,
    tmpdir=None,
    assembler_args=None,
    method="spades",
    memory=None,
    minlen=75,
    screen_accessions=None,
    screen_urls=None,
    mincontiglen=500,
    workdir=None,
    sourdb=None,
    mincovpct=5,
    debug=False,
    pipe=False,
    **kwargs,
):
    """Script runs entire AAFTF pipeline."""
    # script to run entire AAFTF pipeline
    args_dict = {k: v for k, v in locals().items() if k != "kwargs"}
    ram = round(0.75 * get_ram())
    if not memory:
        args_dict["memory"] = str(ram)

    # Helper function to create namespace with required defaults
    def create_namespace(options, required_args=None, **extra_args):
        """Create a Namespace with filtered options and required arguments."""
        namespace_dict = {k: v for (k, v) in args_dict.items() if k in options}
        if required_args:
            namespace_dict.update(required_args)
        namespace_dict.update(extra_args)
        return Namespace(**namespace_dict)

    # Helper function to check step output and handle failures
    def check_step_success(output_file, step_name):
        """Check if step completed successfully."""
        if not check_file(output_file):
            logger.info(f"AAFTF {step_name} failed")
            sys.exit(1)
        return True

    # run trimming with bbduk
    if not check_file(basename + "_1P.fastq.gz"):
        trim_opts = ["memory", "left", "right", "basename", "cpus", "debug", "minlen"]
        trim_args = create_namespace(trim_opts, required_args={"method": "bbduk", "pipe": True, "avgqual": 10})
        trim.run(**vars(trim_args))
    else:
        if right:
            logger.info("AAFTF trim output found: {:} {:}".format(basename + "_1P.fastq.gz", basename + "_2P.fastq.gz"))
        else:
            logger.info("AAFTF trim output found: {:}".format(basename + "_1P.fastq.gz"))
    check_step_success(basename + "_1P.fastq.gz", "trim")

    # run mitochondrial assembly on bbduk trimmed reads
    if right:
        if not check_file(basename + ".mito.fasta"):
            mito_opts = ["left", "right", "out", "minlen", "maxlen", "seed", "starting", "workdir", "pipe", "reference", "memory", "debug"]
            mito_args = create_namespace(
                mito_opts,
                required_args={
                    "left": basename + "_1P.fastq.gz",
                    "right": basename + "_2P.fastq.gz",
                    "out": basename + ".mito.fasta",
                    "minlen": 10000,
                    "maxlen": 100000,
                    "pipe": True,
                    "memory": int(args_dict["memory"]),
                },
                **{x: False for x in mito_opts if x not in args_dict},
            )
            mito.run(**vars(mito_args))
        else:
            logger.info("AAFTF mito output: {}".format(basename + ".mito.fasta"))
    else:
        logger.info("AAFTF mito requires PE reads, " + "skipping mitochondrial de novo assembly")

    # run filtering with bbduk
    if not check_file(basename + "_filtered_1.fastq.gz"):
        filter_opts = ["screen_accessions", "screen_urls", "basename", "cpus", "debug", "memory", "workdir"]
        filter_args = create_namespace(
            filter_opts,
            required_args={
                "aligner": "bbduk",
                "left": basename + "_1P.fastq.gz",
                "right": basename + "_2P.fastq.gz" if right else None,
                "screen_local": [basename + ".mito.fasta"] if check_file(basename + ".mito.fasta") else None,
                "pipe": True,
            },
        )
        aaftf_filter.run(**vars(filter_args))
    else:
        if right:
            logger.info("AAFTF filter output found: {:} {:}".format(basename + "_filtered_1.fastq.gz", basename + "_filtered_2.fastq.gz"))
        else:
            logger.info("AAFTF filter output found: {:}".format(basename + "_filtered_1.fastq.gz"))
    check_step_success(basename + "_filtered_1.fastq.gz", "filter")

    # run assembly with specified method
    assembly_method = method or "spades"
    assembly_file = basename + f".{assembly_method}.fasta"
    if not check_file(assembly_file):
        assemble_opts = ["memory", "cpus", "debug", "workdir", "method", "assembler_args", "tmpdir"]
        asm_extra = {
            "left": basename + "_filtered_1.fastq.gz",
            "right": basename + "_filtered_2.fastq.gz" if right else None,
            "out": assembly_file,
            "pipe": True,
            "method": assembly_method,
            "merged": False,
        }
        if assembly_method == "spades":
            asm_extra["isolate"] = False
            asm_extra["careful"] = True
        asm_args = create_namespace(assemble_opts, required_args=asm_extra)
        assemble.run(**vars(asm_args))
    else:
        logger.info(f"AAFTF assemble output found: {assembly_file}")
    check_step_success(assembly_file, "assemble")

    # run vecscreen
    vecscreen_file = basename + ".vecscreen.fasta"
    if not check_file(vecscreen_file):
        vec_opts = ["cpus", "debug", "workdir"]
        vec_args = create_namespace(vec_opts, required_args={"percent_id": False, "stringency": "high", "infile": assembly_file, "outfile": vecscreen_file, "pipe": True})
        vecscreen.run(**vars(vec_args))
    else:
        logger.info(f"AAFTF vecscreen output found: {vecscreen_file}")
    check_step_success(vecscreen_file, "vecscreen")

    # run sourmash purge
    sourpurge_file = basename + ".sourpurge.fasta"
    if not check_file(sourpurge_file):
        sour_opts = ["cpus", "debug", "workdir", "phylum", "sourdb", "mincovpct"]
        sour_args = create_namespace(
            sour_opts,
            required_args={
                "left": basename + "_filtered_1.fastq.gz",
                "right": basename + "_filtered_2.fastq.gz" if right else None,
                "input": vecscreen_file,
                "outfile": sourpurge_file,
                "kmer": "31",
                "taxonomy": False,
                "pipe": True,
                "sourdb_type": "gbk",
            },
        )
        sourpurge.run(**vars(sour_args))
    else:
        logger.info(f"AAFTF sourpurge output found: {sourpurge_file}")
    check_step_success(sourpurge_file, "sourpurge")

    # run remove duplicates
    rmdup_file = basename + ".rmdup.fasta"
    if not check_file(rmdup_file):
        rmdup_opts = ["cpus", "debug", "workdir"]
        rmdup_args = create_namespace(rmdup_opts, required_args={"input": sourpurge_file, "out": rmdup_file, "minlen": mincontiglen, "percent_id": 95, "percent_cov": 95, "exhaustive": False, "pipe": True})
        rmdup.run(**vars(rmdup_args))
    else:
        logger.info(f"AAFTF rmdup output found: {rmdup_file}")
    check_step_success(rmdup_file, "rmdup")

    # run polish to error-correct
    polish_file = basename + ".polish.fasta"
    if not check_file(polish_file):
        polish_opts = ["cpus", "debug", "workdir", "memory"]
        polish_args = create_namespace(
            polish_opts,
            required_args={
                "method": "pypolca",
                "infile": rmdup_file,
                "outfile": polish_file,
                "left": basename + "_filtered_1.fastq.gz",
                "right": basename + "_filtered_2.fastq.gz" if right else None,
                "longreads": None,
                "pipe": True,
                # pipeline-level --memory is str (matches assemble's spades usage);
                # polish.py divides it by cpus expecting int (see polish.py memperthread).
                "memory": int(args_dict["memory"]),
            },
        )
        polish.run(**vars(polish_args))
    else:
        logger.info(f"AAFTF polish output found: {polish_file}")
    check_step_success(polish_file, "polish")

    # sort and rename
    final_file = basename + ".final.fasta"
    if not check_file(final_file):
        sort_opts = ["debug"]
        sort_args = create_namespace(sort_opts, required_args={"input": polish_file, "out": final_file, "name": "scaffold", "minlen": mincontiglen, "pipe": True})
        aaftf_sort.run(**vars(sort_args))
    else:
        logger.info(f"AAFTF sort output found: {final_file}")
    check_step_success(final_file, "sort")

    # assess the assembly
    assess_opts = ["debug"]
    assess_args = create_namespace(assess_opts, required_args={"input": final_file, "report": False, "telomere_monomer": "TAA[C]+", "telomere_n_repeat": 2, "pipe": True})
    assess.run(**vars(assess_args))

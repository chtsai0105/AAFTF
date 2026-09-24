"""Run the sourmash fast matching kmer tool to look for obvious contaminants."""

import logging
import shutil
from pathlib import Path
from typing import Any

from Bio import SeqIO

from aaftf.utility import align_to_sorted_bam, calc_nx, check_file, cleanup_workdir, execute, fasta_stats, filter_fasta, make_workdir, next_step_name, require_databases, run_cmd

__all__ = ["run"]


logger = logging.getLogger(__name__)


# logging - we may need to think about whether this has
# separate name for the different runfolder
# flake8: noqa: C901
def run(
    input: str,
    outfile: str,
    phylum: list[str],
    workdir: str | None = None,
    cpus: int = 1,
    left: str | None = None,
    right: str | None = None,
    sourdb: str | None = None,
    sourdb_type: str = "gbk",
    kmer: str = "31",
    mincovpct: int = 5,
    taxonomy: bool = False,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run the sourpurge routines to detect and remove contaminant contigs.

    Each contig is classified with ``sourmash lca classify``; classified contigs whose lineage
    contains none of ``phylum`` are dropped. If reads are given they are mapped back with BWA and
    contigs whose mean coverage is at most ``mincovpct`` percent of the average N50-contig
    coverage are dropped too. The sourmash CSV is copied next to ``input`` as
    ``<input stem>.sourmash-taxonomy.csv``.

    Args:
        input: Assembly FASTA to screen.
        outfile: Output FASTA of retained contigs.
        phylum: Taxon names (usually phyla) to keep.
        workdir: Working directory; a temporary one is created when None.
        cpus: Number of threads for BWA.
        left: Forward reads FASTQ; enables the coverage filter.
        right: Reverse reads FASTQ for paired-end data.
        sourdb: Path to a sourmash LCA database; when None the one for ``sourdb_type`` is used.
        sourdb_type: Installed database to use: ``gbk``, ``gtdbrep`` or ``gtdb``.
        kmer: k-mer size (as a string) used for ``sourmash compute``.
        mincovpct: Coverage cutoff as a percent of the average N50-contig coverage.
        taxonomy: Only report the taxonomic classifications and return without filtering.
        debug: Print per-contig details and keep the work directory.
        pipe: Suppress the "next command" hint (set when run from ``pipeline``).
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.
    """
    workdir, custom_workdir = make_workdir(workdir, "sourpurge")
    bamthreads = min(cpus, 4)

    # find reads
    forward_reads, reverse_reads = (None,) * 2
    if left:
        forward_reads = str(Path(left).resolve())
    if right:
        reverse_reads = str(Path(right).resolve())
    if not forward_reads:
        logger.info("Unable to located FASTQ raw reads, low coverage will be skipped. Provide -l,--left (and if paired -r,--right) to enable low coverage filtering.")

    # parse database locations
    if not sourdb:
        # --sourdb_type is restricted to these values by the "sourpurge" menu
        # in _menu.py (choices=["gbk", "gtdbrep", "gtdb"]).
        database = {"gbk": "sm_gbk", "gtdbrep": "sm_gtdbrep", "gtdb": "sm_gtdb"}[sourdb_type]
        sour_db = require_databases([database], hint="or pass --sourdb PATH")[0]
    else:
        sour_db = str(Path(sourdb).resolve())

    # hard coded tmpfile
    assembly_working = "assembly.fasta"
    blob_bam = "remapped.bam"
    shutil.copyfile(input, str(Path(workdir, assembly_working)))
    num_seqs, assembly_size = fasta_stats(str(Path(workdir, assembly_working)))
    logger.info(f"Assembly is {num_seqs:,} contigs and {assembly_size:,} bp")

    # now filter for taxonomy with sourmash lca classify
    logger.info("Running SourMash to get taxonomy classification for each contig")
    sour_sketch = Path(assembly_working).name + ".sig"

    sour_compute = ["sourmash", "compute", "-k", kmer, "--scaled=1000", "--singleton", assembly_working]
    run_cmd(sour_compute, debug, cwd=workdir)
    sour_classify = ["sourmash", "lca", "classify", "--db", sour_db, "--query", sour_sketch]
    # output csv: ID,status,superkingdom,phylum,class,order,family,genus,species,strain
    contig_taxonomy = {}
    unique_tax: set[str] = set()
    sourmash_tsv = str(Path(workdir, "sourmash.csv"))
    with open(sourmash_tsv, "w") as sour_out:
        for line in execute(sour_classify, workdir, debug):
            sour_out.write(line)
            if not line or line.startswith("\n") or line.startswith("ID") or line.count(",") < 9:
                continue
            line = line.strip()
            cols = line.split(",")
            if "found" in cols[1]:
                idx = 1
                contig_taxonomy[cols[0]] = cols[idx + 1 :]
                tax_clean = [x for x in cols[idx + 1 :] if x]
                unique_tax.add(";".join(tax_clean))
            elif cols[1].strip() == "nomatch":
                idx = cols.index("nomatch")
                contig_taxonomy[cols[0]] = cols[idx + 1 :]
    logger.info("Found {:} taxonomic classifications for contigs:\n{:}".format(len(unique_tax), "\n".join(unique_tax)))
    if taxonomy:  # --taxonomy only reports the classifications
        return
    tax_to_drop = []
    for k, v in contig_taxonomy.items():
        v = [x for x in v if x]  # remove empty items from list
        if debug:
            logger.info(f"{k}\t{v}")
        if len(v) > 0:
            if not any(i in v for i in phylum):
                tax_to_drop.append(k)

    # drop contigs from taxonomy before calculating coverage
    logger.info(f"Dropping {len(tax_to_drop)} contigs from taxonomy screen")
    sour_tax = str(Path(workdir, "sourmashed-tax-screen.fasta"))
    tax_drop = set(tax_to_drop)
    filter_fasta(str(Path(workdir, assembly_working)), sour_tax, lambda seq_id: seq_id not in tax_drop)

    # only do coverage trimming if reads provided
    contigs_to_drop = []  # this will be empty if no reads given to gather by coverage
    if forward_reads:
        # check if BAM present, if so skip running
        if not Path(workdir, blob_bam).is_file():
            # index
            bwa_index = ["bwa", "index", Path(sour_tax).name]
            logger.info("Building BWA index")
            run_cmd(bwa_index, debug, cwd=workdir)
            # mapped reads to assembly using BWA
            bwa_cmd = [
                "bwa",
                "mem",
                "-t",
                str(cpus),
                Path(sour_tax).name,  # assembly index base
                forward_reads,
            ]
            if reverse_reads:
                bwa_cmd.append(reverse_reads)

            logger.info("Aligning reads to assembly with BWA")
            align_to_sorted_bam(bwa_cmd, str(Path(workdir, blob_bam)), bamthreads, cwd=workdir, debug=debug)

        # now calculate coverage from BAM file
        logger.info("Calculating read coverage per contig")
        fasta_bed = str(Path(workdir, "assembly.bed"))
        lengths = []
        with open(fasta_bed, "w") as bedout:
            with open(sour_tax) as seq_in:
                for record in SeqIO.parse(seq_in, "fasta"):
                    bedout.write(f"{record.id}\t{0}\t{len(record.seq)}\n")
                    lengths.append(len(record.seq))

        n50, _ = calc_nx(lengths)
        coverage = {}
        coverage_bed = str(Path(workdir, "coverage.bed"))
        cov_cmd = ["samtools", "bedcov", Path(fasta_bed).name, blob_bam]
        with open(coverage_bed, "w") as bed_out:
            for line in execute(cov_cmd, workdir, debug):
                bed_out.write(line)

                if not line or line.startswith("\n") or line.count("\t") < 3:
                    continue

                line = line.strip()
                cols = line.split("\t")
                cov = int(cols[3]) / float(cols[2])
                coverage[cols[0]] = (int(cols[2]), cov)

        # get average coverage of N50 contigs
        n50_cov = []
        for contig, (contig_len, cov) in coverage.items():
            if debug:
                print(f"{contig}; Len: {contig_len}; Cov: {cov:.2f}")
            if contig_len >= n50:
                n50_cov.append(cov)
        if not n50_cov:
            logger.warning(f"No contig of at least the N50 length ({n50:,} bp) has coverage data; skipping the low-coverage filter")
        else:
            n50_avg_cov = sum(n50_cov) / len(n50_cov)
            minpct = mincovpct / 100
            # should we make this a variable? 5% was something arbitrary
            min_coverage = float(n50_avg_cov * minpct)
            logger.info(f"Average coverage for N50 contigs is {int(n50_avg_cov)}X")

            # Start list of contigs to drop
            for contig, (_, cov) in coverage.items():
                if cov <= min_coverage:
                    contigs_to_drop.append(contig)
            logger.info(f"Found {len(contigs_to_drop):,} contigs with coverage less than {min_coverage:.2f}X ({mincovpct}%)")

    if debug:
        print("Contigs dropped due to coverage: {:}".format(",".join(contigs_to_drop)))
        print("Contigs dropped due to taxonomy: {:}".format(",".join(tax_to_drop)))

    drop_final = set(contigs_to_drop + tax_to_drop)
    logger.info(f"Dropping {len(drop_final):,} total contigs based on taxonomy and coverage")
    num_seqs, assembly_size = filter_fasta(sour_tax, outfile, lambda seq_id: seq_id not in drop_final)
    logger.info(f"Sourpurged assembly is {num_seqs:,} contigs and {assembly_size:,} bp")
    next_out = next_step_name(outfile, ".rmdup.fasta")

    if check_file(sourmash_tsv):
        baseinput = Path(input).name
        basedir = str(Path(input).parent)
        if "." in baseinput:
            baseinput = baseinput.rsplit(".", 1)[0]

        shutil.copy(sourmash_tsv, str(Path(basedir, baseinput + ".sourmash-taxonomy.csv")))

    cleanup_workdir(workdir, debug, custom_workdir)

    if not pipe:
        logger.info(f"Your next command might be:\nAAFTF rmdup -i {outfile} -o {next_out}")

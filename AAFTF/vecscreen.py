"""This runs routines to identify contaminant contigs.

The contaminants are presumably sequences that were not screened out
in the filter step.
The vector library UniVec and known or user specified contaminanting
sequences (by GenBank accession number) can be provided for
additional cleanup

The default libraries for screening are located in resources.py
and include common Euk, Prok, and MITO contaminants.
"""

import csv
import logging
import os
from pathlib import Path
from subprocess import DEVNULL, call

# biopython needed
from Bio import SeqIO

from AAFTF.utility import cleanup_workdir, concat_files, make_workdir, next_step_name, printCMD, require_databases, write_fasta

logger = logging.getLogger(__name__)

BlastPercent_ID_ContamMatch = "90.0"
BlastPercent_ID_MitoMatch = "98.6"

# VecScreen matches
"""
Strong Match to Vector
(Expect 1 random match in 1,000,000 queries of length 350 kb.)
Terminal match with Score ≥ 24.
Internal match with Score ≥ 30.
Moderate Match to Vector
(Expect 1 random match in 1,000 queries of length 350 kb.)
Terminal match with Score 19 to 23.
Internal match with Score 25 to 29.
Weak Match to Vector
(Expect 1 random match in 40 queries of length 350 kb.)
Terminal match with Score 16 to 18.
Internal match with Score 23 to 24.
Segment of Suspect Origin
Any segment of fewer than 50 bases between two vector matches
    or between a match and an end.
"""


def run(
    infile,
    outfile,
    workdir=None,
    cpus=1,
    percent_id=None,
    stringency="high",
    debug=False,
    pipe=False,
    **kwargs,
):
    """Runs vectorscreening via BLASTN against a vectorDB.

    Pipeline: build the contamination BLAST databases, screen out Euk/Prok
    contaminants, screen for mitochondrial contigs, run repeated VecScreen
    (UniVec) rounds to trim/split vector hits, then write the cleaned
    assembly and a separate mitochondrial-contigs FASTA.
    """
    workdir, custom_workdir = make_workdir(workdir, "vecscreen")
    percentid_cutoff = percent_id or BlastPercent_ID_ContamMatch

    # final_outfile/outdir/prefix are derived from the user's --outfile once,
    # up front, so nothing later in the pipeline can accidentally clobber them
    # (each stage below only sees its own function-local variables).
    final_outfile = Path(outfile).name
    outdir = str(Path(final_outfile).parent)
    if ".f" in final_outfile:
        prefix = final_outfile.rsplit(".f", 1)[0]
    else:
        prefix = str(os.getpid())
    if not final_outfile:
        final_outfile = f"{prefix}.vecscreen.fasta"
    _build_contam_databases(workdir)

    contigs_to_remove = {}
    regions_to_trim = _screen_euk_prok_contamination(infile, workdir, prefix, cpus, percentid_cutoff)
    eukCleaned = _write_euk_cleaned(infile, regions_to_trim, workdir, prefix)
    mitoHits = _screen_mitochondria(eukCleaned, workdir, prefix, cpus, contigs_to_remove)
    outfile_vec = _run_vecscreen_rounds(eukCleaned, workdir, prefix, cpus, stringency, contigs_to_remove)

    logger.info(f"{len(contigs_to_remove):,} contigs will be removed:")
    for k, v in sorted(contigs_to_remove.items()):
        print(f"\t{k} --> dbhit={v[0]}; hit={v[1]}; pident={v[2]}")

    # this could instead use the outfile and strip
    # .fasta/fsa/fna and add mito on it I suppose, but assumes
    # a bit about the naming structure
    mitochondria = str(Path(outdir, prefix + ".mitochondria.fasta"))
    _write_final_outputs(outfile_vec, contigs_to_remove, mitoHits, final_outfile, mitochondria)

    nextOut = next_step_name(final_outfile, ".sourpurge.fasta")
    if not pipe:
        logger.info("Your next command might be:\n" + "AAFTF sourpurge -i {:} -o {:} -c {:} --phylum {:}".format(final_outfile, nextOut, cpus, "Ascomycota"))

    cleanup_workdir(workdir, debug, custom_workdir)


# BLAST database name used by the screens -> `AAFTF download` database it is built from
_CONTAM_BLAST_DBS = {"UniVec": "univec", "CONTAM_EUKS": "euks", "CONTAM_PROKS": "proks", "MITO": "mitodb"}


def _build_contam_databases(workdir):
    """Build a BLAST nucleotide DB in ``workdir`` for each contamination screen, from the downloaded databases."""
    logger.info("Building BLAST databases for contamination screen.")
    sources = require_databases(list(_CONTAM_BLAST_DBS.values()))
    for blast_name, source in zip(_CONTAM_BLAST_DBS, sources):
        combined_fasta = str(Path(workdir, f"{blast_name}.fasta"))
        concat_files([source], combined_fasta)
        _make_blastdb("nucl", combined_fasta, str(Path(workdir, blast_name)))


def _make_blastdb(type, file, name):
    """Create the BLASTN database for the vecscreen vector search."""
    idxfile = name
    if type == "nucl":
        idxfile += ".nin"
    else:
        idxfile += ".pin"
    idxexists = Path(idxfile).exists()
    if not idxexists or Path(idxfile).stat().st_ctime < Path(file).stat().st_ctime:
        cmd = ["makeblastdb", "-dbtype", type, "-in", file, "-out", name]
        printCMD(cmd)
        call(cmd, stdout=DEVNULL, stderr=DEVNULL)


def _run_blastn_screen(query, workdir, prefix, dbname, cpus, percent_identity):
    """Run blastn (tab-6 output) for ``query`` against a DB built in ``workdir``, and return the parsed rows."""
    blastreport = str(Path(workdir, f"{dbname}.{prefix}.blastn"))
    blastnargs = ["blastn", "-query", query, "-db", str(Path(workdir, dbname)), "-num_threads", str(cpus), "-dust", "yes", "-soft_masking", "true", "-perc_identity", percent_identity, "-lcase_masking", "-outfmt", "6", "-out", blastreport]
    printCMD(blastnargs)
    call(blastnargs)
    with open(blastreport) as report:
        return list(csv.reader(report, delimiter="\t"))


def _screen_euk_prok_contamination(infile, workdir, prefix, cpus, percentid_cutoff):
    """BLASTN infile against the Euk/Prok contamination DBs.

    Returns {contig_id: [(start, end, db_name, hit_id, pident), ...]}.
    """
    # qaccver saccver pident length mismatch gapopen qstart qend
    # sstart send evalue bitscore
    regions_to_trim = {}
    for contam in ["CONTAM_EUKS", "CONTAM_PROKS"]:
        logger.info(f"{contam} Contamination Screen")
        for row in _run_blastn_screen(infile, workdir, prefix, contam, cpus, percentid_cutoff):
            if (float(row[2]) >= 98.0 and int(row[3]) >= 50) or (float(row[2]) >= 94.0 and int(row[3]) >= 100) or (float(row[2]) >= 90.0 and int(row[3]) >= 200):
                start, end = sorted([int(row[6]), int(row[7])])
                regions_to_trim.setdefault(row[0], []).append((start, end, contam, row[1], float(row[2])))
        logger.info(f"{contam} screening finished")
    return regions_to_trim


def _write_euk_cleaned(infile, regions_to_trim, workdir, prefix):
    """Split out Euk/Prok-contaminated regions found by ``_screen_euk_prok_contamination``.

    Writes a cleaned FASTA and returns its path, or returns ``infile``
    unchanged if there was nothing to trim.
    """
    if not regions_to_trim:
        return infile

    eukCleaned = str(Path(workdir, f"{prefix}.euk-prot_cleaned.fasta"))
    with open(eukCleaned, "w") as cleanout, open(infile) as fastain:
        for record in SeqIO.parse(fastain, "fasta"):
            if record.id not in regions_to_trim:
                write_fasta(cleanout, record.id, str(record.seq))
            else:
                Seq = str(record.seq)
                regions = regions_to_trim[record.id]
                logger.info(f"Splitting {record.id} for contamination: {regions}")
                lastpos = 0
                for i, x in enumerate(regions):
                    # x[0] is a 1-based BLAST start; subtract one so the
                    # slice end doesn't retain the first contaminant base.
                    newSeq = Seq[lastpos : x[0] - 1]
                    lastpos = x[1]
                    write_fasta(cleanout, f"split{i}_{record.id}", newSeq)
                    if i == len(regions) - 1:
                        newSeq = Seq[x[1] :]
                        write_fasta(cleanout, f"split{i + 1}_{record.id}", newSeq)
    return eukCleaned


def _screen_mitochondria(eukCleaned, workdir, prefix, cpus, contigs_to_remove):
    """BLASTN against the MITO DB; flags long hits for removal in ``contigs_to_remove`` (mutated in place).

    Returns the list of contig ids identified as mitochondrial.
    """
    logger.info("Mitochondria Contamination Screen")
    mitoHits = []
    for row in _run_blastn_screen(eukCleaned, workdir, prefix, "MITO", cpus, BlastPercent_ID_MitoMatch):
        if int(row[3]) >= 120:
            contigs_to_remove[row[0]] = ("MitoScreen", row[1], float(row[2]))
            mitoHits.append(row[0])
    logger.info("Mito screening finished.")
    return mitoHits


def _run_vecscreen_rounds(eukCleaned, workdir, prefix, cpus, stringency, contigs_to_remove):
    """Repeatedly BLASTN against UniVec, trimming/splitting vector hits each round until none remain.

    ``contigs_to_remove`` is mutated in place by ``_parse_clean_blastn``.
    Returns the path to the final, fully vector-cleaned FASTA.
    """
    logger.info("Starting VecScreen, will remove terminal matches and split internal matches")
    rnd = 0
    count = 1
    cleanfile = eukCleaned
    while count > 0:
        filepref = f"{prefix}.r{rnd}"
        report = str(Path(workdir, f"{filepref}.vecscreen.tab"))
        if not Path(report).exists():
            cmd = [
                "blastn",
                "-task",
                "blastn",
                "-reward",
                "1",
                "-penalty",
                "-5",
                "-gapopen",
                "3",
                "-gapextend",
                "3",
                "-dust",
                "yes",
                "-soft_masking",
                "true",
                "-evalue",
                "700",
                "-searchsp",
                "1750000000000",
                "-db",
                str(Path(workdir, "UniVec")),
                "-outfmt",
                "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore score qlen",
                "-num_threads",
                str(cpus),
                "-query",
                eukCleaned,
                "-out",
                report,
            ]
            printCMD(cmd)
            call(cmd)
        logger.info(f"Parsing VecScreen round {rnd + 1}: {filepref} for {report}")

        (count, cleanfile) = _parse_clean_blastn(eukCleaned, str(Path(workdir, filepref)), report, stringency, contigs_to_remove)
        logger.info(f"count is {count} cleanfile is {cleanfile}")
        if count > 0:  # vector matches remain; re-screen the newly trimmed/split sequences
            rnd += 1
            eukCleaned = cleanfile

    return cleanfile


def _write_final_outputs(outfile_vec, contigs_to_remove, mitoHits, outfile, mitochondria):
    """Split the vecscreen output into the cleaned assembly and a mitochondrial-contigs FASTA."""
    n_clean = n_mito = 0
    with open(outfile, "w") as oh, open(mitochondria, "w") as mh:
        for record in SeqIO.parse(outfile_vec, "fasta"):
            if record.id not in contigs_to_remove:
                write_fasta(oh, record.description, str(record.seq))
                n_clean += 1
            elif record.id in mitoHits:
                write_fasta(mh, record.description, str(record.seq))
                n_mito += 1
    logger.info(f"Writing {n_clean:,} cleaned contigs to: {outfile}")
    logger.info(f"Writing {n_mito:,} mitochondrial contigs to: {mitochondria}")


def _parse_clean_blastn(fastafile, prefix, blastn, stringent, contigs_to_remove=None):
    """Parse a VecScreen BLASTN report and write a vector-trimmed/split FASTA.

    Returns (found_vector_seq, cleaned_fasta_path).
    """
    cleaned = prefix + ".clean.fsa"
    if contigs_to_remove is None:
        contigs_to_remove = {}
    vec_hits, found_vector_seq = _classify_vector_hits(blastn, stringent, contigs_to_remove)
    _write_trimmed_and_split(fastafile, vec_hits, cleaned)
    return found_vector_seq, cleaned


def _classify_vector_hits(blastn, stringent, contigs_to_remove):
    """Read a VecScreen-style BLASTN tab report and classify each surviving hit.

    Input rows have columns: qaccver saccver pident length mismatch gapopen
    qstart qend sstart send evalue bitscore score qlen

    Each hit is classified as terminal (within 25bp of either end of the
    query) or internal, then scored weak/moderate/strong per the VecScreen
    table (see module docstring):
    https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/#Moderate
    (using ``score``, not ``bitscore``). Hits on a contig already in
    ``contigs_to_remove`` are skipped. Weak hits are always dropped;
    ``stringent == "high"`` additionally keeps moderate hits, otherwise only
    strong hits are kept.

    Returns (vec_hits, found_vector_seq):
      vec_hits: {contig_id: [(hit_id, qlen, loc, score, terminal, position), ...]}
      found_vector_seq: number of qualifying hits found
    """
    vec_hits = {}
    found_vector_seq = 0
    with open(blastn) as vectab:
        rdr = csv.reader(vectab, delimiter="\t")
        for row in rdr:
            (qaccver, saccver, pid, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, score, qlen) = row
            if qaccver in contigs_to_remove:
                continue

            loc = sorted([int(qstart), int(qend)])
            terminal = False
            position = None
            if loc[0] <= 25:
                terminal = True
                position = "5"
            if (int(qlen) - loc[1]) <= 25:
                terminal = True
                position = "3"

            score = int(score)
            match_strength = 0  # weak=0, moderate=1, strong=2
            if terminal:
                if score >= 19:
                    match_strength = 1
                if score >= 24:
                    match_strength = 2
            else:
                if score >= 25:
                    match_strength = 1
                if score >= 30:
                    match_strength = 2

            if match_strength == 0:
                continue
            if stringent != "high" and match_strength < 2:
                continue

            found_vector_seq += 1
            vec_hits.setdefault(qaccver, []).append((saccver, int(qlen), loc, score, terminal, position))
    return vec_hits, found_vector_seq


# flake8: noqa: C901
def _write_trimmed_and_split(fastafile, vec_hits, cleaned):
    """Write ``cleaned`` FASTA, trimming terminal vector hits and splitting out internal ones.

    For each record with hits in ``vec_hits``: terminal 5' hits push the
    kept region's start forward, terminal 3' hits pull its end back, and
    internal hits carve the contig into multiple ``splitN_<id>`` pieces
    around the vector regions (any hits within 50bp of each other collapse
    into one cut, handled naturally across ``_run_vecscreen_rounds``'
    repeated rounds rather than here). Records/pieces shorter than 200bp are
    dropped.
    """
    with open(cleaned, "w") as output_handle:
        for record in SeqIO.parse(fastafile, "fasta"):
            seq_str = str(record.seq)
            if record.id not in vec_hits:
                if len(record.seq) >= 200:
                    write_fasta(output_handle, record.id, seq_str)
                continue

            five_end = 0
            three_end = len(record.seq)
            internal_regions = []
            for hit_id, qlen, loc, score, terminal, position in vec_hits[record.id]:
                if terminal and position == "5":
                    five_end = max(five_end, loc[1])
                elif terminal and position == "3":
                    # loc[0] (qstart) is 1-based; the slice end must be
                    # one less so the first vector/contaminant base isn't kept.
                    three_end = min(three_end, loc[0] - 1)
                elif loc not in internal_regions:
                    internal_regions.append(loc)

            sorted_internals = sorted(internal_regions, key=lambda region: region[0])
            slicer = [five_end]
            for start, end in sorted_internals:
                # start (qstart) is 1-based and used as a slice end below, so
                # shift by one to avoid retaining the first hit base. end
                # (qend) is used as the following slice start, already correct.
                slicer += [start - 1, end]
            slicer.append(three_end)
            keep_regions = list(_group(slicer, 2))

            if len(keep_regions) < 2:
                logger.info(f"Terminal trimming {record.id} to {keep_regions}")
                newSeq = seq_str[keep_regions[0][0] : keep_regions[0][1]]
                if len(newSeq) >= 200:
                    write_fasta(output_handle, record.id, newSeq)
            else:
                logger.info(f"Splitting contig {record.id} into {keep_regions}")
                for num, (start, end) in enumerate(keep_regions):
                    newSeq = seq_str[start:end]
                    if len(newSeq) >= 200:
                        write_fasta(output_handle, f"split{num + 1}_{record.id}", newSeq)


def _group(lst, n):
    """This groups sets by a size."""
    for i in range(0, len(lst), n):
        val = lst[i : i + n]
        if len(val) == n:
            yield tuple(val)

"""URLS and hardcoded resource links."""

__all__ = ["NCBI", "CONTAMINANT_ACCESSIONS", "DB_LINKS", "EUTILS", "SEQ_DBS", "MITO_SEQS", "FCSADAPTOR", "DATABASES"]


NCBI = "https://ftp.ncbi.nlm.nih.gov"

CONTAMINANT_ACCESSIONS = {"phiX": [f"{NCBI}/genomes/all/GCF/000/819/615/" + "GCF_000819615.1_ViralProj14015/" + "GCF_000819615.1_ViralProj14015_genomic.fna.gz"]}

DB_LINKS = {
    "UniVec": [f"{NCBI}/pub/UniVec/UniVec"],
    "CONTAM_EUKS": [f"{NCBI}/pub/kitts/contam_in_euks.fa.gz"],
    "CONTAM_PROKS": [f"{NCBI}/pub/kitts/contam_in_prok.fa"],
    "MITO": [f"{NCBI}/refseq/release/mitochondrion/" + "mitochondrion.1.1.genomic.fna.gz"],
    "sourmash_gbk": [{"version": "2017.11.07", "filename": "genbank-k31.lca.json.gz", "url": "https://osf.io/4f8n3/download"}],
    # first in list is default, will fix someday to allow choosing the version
    "sourmash_gtdbrep": [
        {"version": "rs220", "filename": "gtdb-rs220-reps.k31.lca.json.gz", "url": "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs220/gtdb-reps-rs220-k31.dna.zip"},  # noqa: E501
    ],
    # first in list is default, will fix someday to allow choosing the version
    "sourmash_gtdb": [
        {"version": "rs220", "filename": "gtdb-rs220-k31.lca.json.gz", "url": "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs220/gtdb-rs220-k31.dna.zip"},  # noqa: E501
    ],
}

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# URL templates for fetching one sequence by accession (``% accession``)
SEQ_DBS = {
    "nucleotide": f"{EUTILS}/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta",
    "nucleotide_ebi": "https://www.ebi.ac.uk/ena/data/view/%s?display=fasta",
    "nucleotide_ncbi": f"{EUTILS}/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta",
}

MITO_SEQS = {
    "COB1": (
        "atgagaattttaaaaagtcatcctttattaaaattagttaatagttatattattg"
        + "attcaccacaaccttctaatattagttatttatgaaattttggatctttattagc"
        + "tttatgtttagttatacaaattgtaactggtgttacattagctatgcactataca"
        + "cctaatgttgatttagcttttaattctgtagaacatattatgagagatgtaaata"
        + "atggttgattaataagatatttacatgctaatactgcttcagcattctttttctt"
        + "agttatatttacatataggtagaggattatattatggttcatataaatcacctag"
        + "aacttaacatgagctattgg"
    ),
}

"""NCBI Foreign Contaminant Screen tool links"""

FCSADAPTOR = {
    "VERSION": "0.5.5",
    "SIF": "fcs-adaptor.sif",
    "SIFLOCAL": "fcs-adaptor.%s.sif",
    "SIFURL": "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/",
    "EXEURL": "https://raw.githubusercontent.com/ncbi/fcs/v%s/dist/run_fcsadaptor.sh",  # noqa: E501
    "DOCKERIMAGE": "ncbi/fcs-adaptor:%s",
}

# Databases `AAFTF download` can fetch, keyed by abbreviation. `used_by` lists the subcommands that
# read each one; `executable` marks files that must be runnable after download.
DATABASES = {
    "phix": {"filename": CONTAMINANT_ACCESSIONS["phiX"][0].rsplit("/", 1)[-1], "url": CONTAMINANT_ACCESSIONS["phiX"][0], "used_by": "filter"},
    "univec": {"filename": DB_LINKS["UniVec"][0].rsplit("/", 1)[-1], "url": DB_LINKS["UniVec"][0], "used_by": "filter, vecscreen"},
    "euks": {"filename": DB_LINKS["CONTAM_EUKS"][0].rsplit("/", 1)[-1], "url": DB_LINKS["CONTAM_EUKS"][0], "used_by": "vecscreen"},
    "proks": {"filename": DB_LINKS["CONTAM_PROKS"][0].rsplit("/", 1)[-1], "url": DB_LINKS["CONTAM_PROKS"][0], "used_by": "vecscreen"},
    "mitodb": {"filename": DB_LINKS["MITO"][0].rsplit("/", 1)[-1], "url": DB_LINKS["MITO"][0], "used_by": "vecscreen"},
    "sm_gbk": {**{k: DB_LINKS["sourmash_gbk"][0][k] for k in ("filename", "url")}, "used_by": "sourpurge (--sourdb_type gbk)"},
    "sm_gtdbrep": {**{k: DB_LINKS["sourmash_gtdbrep"][0][k] for k in ("filename", "url")}, "used_by": "sourpurge (--sourdb_type gtdbrep)"},
    "sm_gtdb": {**{k: DB_LINKS["sourmash_gtdb"][0][k] for k in ("filename", "url")}, "used_by": "sourpurge (--sourdb_type gtdb)"},
    "fcs_script": {"filename": "run_fcsadaptor.sh", "url": FCSADAPTOR["EXEURL"] % FCSADAPTOR["VERSION"], "used_by": "fcs_screen", "executable": True},
    "fcs_image": {
        "filename": FCSADAPTOR["SIFLOCAL"] % FCSADAPTOR["VERSION"],
        "url": "/".join([FCSADAPTOR["SIFURL"].rstrip("/"), FCSADAPTOR["VERSION"], FCSADAPTOR["SIF"]]),
        "used_by": "fcs_screen (singularity)",
    },
}

============
Installation
============

AAFTF can be run three ways: a **conda/pip environment** you assemble yourself, a **pixi**-managed
environment (the same mechanism used to build the project's own Docker/Singularity images), or
directly from a pre-built **Docker/Singularity/Apptainer container**. All external bioinformatics
tools (SPAdes, BBTools, BLAST+, sourmash, bwa, minimap2, pypolca, polypolish, NextPolish2, mosdepth, ...)
must be reachable on ``$PATH`` -- the three methods below differ only in *how* that PATH gets populated.

Option 1: conda environment
============================

.. code-block:: bash

    conda create -n aaftf -c bioconda -c conda-forge "python>=3.10,<3.15" \
        bbmap trimmomatic bowtie2 bwa freebayes yak sourmash blast minimap2 \
        spades megahit novoplasty "biopython>=1.88" fastp pypolca polypolish \
        nextpolish2 racon unicycler mosdepth "matplotlib>=3" "samtools>=1.24"

    conda activate aaftf
    pip install AAFTF
    # or track the latest from GitHub:
    python -m pip install git+https://github.com/stajichlab/AAFTF.git

.. warning::
   AAFTF requires samtools >= 1.0 (ideally >= 1.24 for best performance -- newer samtools
   avoids writing unsorted temp BAM files to disk). If you assemble your own environment,
   watch for bioconda packages that bundle an old, conflicting samtools binary of their own
   (older MaSuRCA builds were a historical example of this) -- AAFTF's own dependency lists
   (``pyproject.toml``, ``environment.yml``) no longer pull in any such package.

Option 2: pixi (recommended for reproducible/locked environments)
====================================================================

The repository ships a ``pyproject.toml`` (``[tool.pixi.*]`` tables) / ``pixi.lock`` pinning every dependency (this is what the
Docker and Singularity images are built from). From a checkout of the repository:

.. code-block:: bash

    curl -fsSL https://pixi.sh/install.sh | bash   # install pixi, if needed
    cd AAFTF
    pixi install                       # creates the "default" (editable, dev) environment
    pixi run install-bowtie2           # rebuild bowtie2 from source (fixes an AVX2/x86-64-v3
                                        # runtime fallback bug in the conda binary)
    pixi shell                         # activate the environment
    AAFTF --version

To install a specific tagged release instead of an editable checkout, use the ``release``
environment (tracks the ``aaftf`` package pinned in ``pyproject.toml`` (``[tool.pixi.*]`` tables)):

.. code-block:: bash

    pixi install --environment release
    pixi run --environment release AAFTF --version

Option 3: Docker / Singularity / Apptainer
============================================

Prebuilt containers bundle every dependency, including ``pypolca`` (POLCA-style polishing that
works with modern samtools) and a source-built ``bowtie2``. Only ``fcs_screen`` additionally requires a
*separate* container engine (singularity/apptainer or docker) at runtime, since NCBI FCS-adaptor
itself ships as a container image invoked from inside AAFTF -- see :doc:`commands/fcs_screen`.

**Docker**

.. code-block:: bash

    # Build (AAFTF_VERSION must be supplied explicitly; .git is excluded from the build context)
    docker build --build-arg AAFTF_VERSION=$(git describe --tags --always) -t aaftf:latest .

    # Run
    docker run --rm aaftf:latest AAFTF --help

    # Run with data + a persistent reference-database directory mounted in
    docker run --rm \
        -v /path/to/data:/data \
        -v /path/to/aaftf_db:/opt/aaftf_db \
        aaftf:latest AAFTF trim --left /data/R1.fq.gz --right /data/R2.fq.gz

A second, simpler conda-based image is available via ``Dockerfile.conda`` for environments where
the pixi-based build is undesirable; build/run commands are the same.

**Singularity / Apptainer**

.. code-block:: bash

    # Build (defaults: git_ref=main, skip_db_download=1 -- skips baking the multi-GB sourmash DB
    # into the image; bind-mount or download AAFTF_DB separately at run time instead)
    singularity build AAFTF.sif AAFTF.def

    # Build a specific tagged release
    singularity build --build-arg git_ref=v0.7.0 --build-arg skip_db_download=1 AAFTF.sif AAFTF.def

    # Run
    singularity exec AAFTF.sif AAFTF --help
    singularity exec --bind /path/to/data:/data AAFTF.sif \
        AAFTF trim --left /data/R1.fq.gz --right /data/R2.fq.gz

    # Or use the built-in runscript
    singularity run --bind /path/to/data:/data AAFTF.sif trim --left /data/R1.fq.gz --right /data/R2.fq.gz

The database directory defaults to ``/opt/aaftf_db`` inside the container (override with
``AAFTF_DB``). Both the Docker entrypoint and the Singularity ``%environment``/``/etc/profile.d``
hook ensure the pixi-managed tool PATH survives both interactive (``singularity exec``) and login
shells (``bash -l``, as used by SLURM/Nextflow task scripts).

Setting up the reference database (``AAFTF_DB``)
=====================================================

Several subcommands (``filter``, ``vecscreen``, ``sourpurge``, ``fcs_screen``) need UniVec, PhiX,
eukaryotic/prokaryotic/mitochondrial contamination screens, sourmash taxonomy databases, and/or
the NCBI FCS-adaptor resources. Download them once with ``AAFTF database``; those subcommands then
read them from the database folder (and tell you what to download if something is missing). Without
``$AAFTF_DB`` this is ``~/.cache/aaftf``; because some databases are several GB, it is best to
point ``$AAFTF_DB`` at a folder with plenty of space (it can also list several folders separated
by ``:``, like ``$PATH``):

.. code-block:: bash

    export AAFTF_DB=/path/with/space/aaftf_db
    AAFTF database                                        # list databases
    AAFTF database phix univec euks proks mitodb sm_gbk   # download some

See :doc:`commands/database` for the list of databases and which subcommands use each (the
sourmash indices are large, so fetch only the one you need).

Runs launched inside the pixi-managed Docker/Singularity images default ``AAFTF_DB`` to
``/opt/aaftf_db``; bind-mount or ``-e AAFTF_DB=...`` to point at your own copy instead of
re-downloading it on every run.

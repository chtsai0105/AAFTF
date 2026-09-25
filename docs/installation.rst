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

    # only what the default pipeline needs (same as the pixi "default" environment);
    # installs AAFTF from PyPI
    conda env create -f environment.yml
    conda activate aaftf

    # or everything AAFTF can use (same as the pixi "complete" environment), with AAFTF
    # installed from this checkout in editable mode
    conda env create -f environment.dev.yml
    conda activate aaftf-dev

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
    pixi install -e complete               # every tool AAFTF can use
    pixi shell -e complete                 # activate the environment
    AAFTF --version

Three environments are defined, each installing AAFTF from the checkout in editable mode:

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Environment
     - Features
     - Contents
   * - ``default``
     - default
     - Only what ``AAFTF pipeline`` needs with its default settings (BBTools, SPAdes, BLAST+,
       sourmash, bwa, samtools, minimap2). ``pixi install`` alone creates this one.
   * - ``complete``
     - default, optional
     - Plus every tool used by optional steps (``polish``, ``depth``, ``mito``, ``fcs_screen``,
       ``fcs_gx_purge``) and non-default options (fastp, Trimmomatic, MEGAHIT, Unicycler, ...).
       The Docker and Singularity images use this environment.
   * - ``dev``
     - default, optional, dev
     - Plus pytest, pytest-cov and pre-commit, for working on AAFTF itself.

``AAFTF dependency`` reports what an environment is missing.

The ``complete`` and ``dev`` environments include bowtie2 (used by ``filter --aligner bowtie2``) from
conda. That build cannot use AVX2 instructions (it falls back from its x86-64-v3 version at
runtime), so it runs slower than it could. If you use bowtie2 and want AVX2 support, rebuild it from
source into the environment once (this replaces the conda binaries; the Docker and Singularity
images already do this):

.. code-block:: bash

    pixi run -e complete install-bowtie2   # or: pixi run -e dev install-bowtie2

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
        aaftf:latest AAFTF trim --read1 /data/R1.fq.gz --read2 /data/R2.fq.gz

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
        AAFTF trim --read1 /data/R1.fq.gz --read2 /data/R2.fq.gz

    # Or use the built-in runscript
    singularity run --bind /path/to/data:/data AAFTF.sif trim --read1 /data/R1.fq.gz --read2 /data/R2.fq.gz

The database directory defaults to ``/opt/aaftf_db`` inside the container (override with
``AAFTF_DB``). Both the Docker entrypoint and the Singularity ``%environment``/``/etc/profile.d``
hook ensure the pixi-managed tool PATH survives both interactive (``singularity exec``) and login
shells (``bash -l``, as used by SLURM/Nextflow task scripts).

Checking the installation
=========================

``AAFTF dependency`` lists every external tool and Python package AAFTF uses, in two groups: those
``AAFTF pipeline`` needs with its default settings (BBTools and Java, SPAdes, BLAST+, sourmash,
bwa, samtools, minimap2; biopython, psutil), and those only needed by optional steps
(``polish``, ``depth``, ``mito``, ``fcs_screen``, ``fcs_gx_purge``) or non-default options
(e.g. ``filter --aligner bowtie2``, ``assemble --method megahit``). It exits with an error only
when something from the first group is missing; ``-q`` prints just the OK/MISSING status.

.. code-block:: bash

    AAFTF dependency

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

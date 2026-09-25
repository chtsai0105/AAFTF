========
pipeline
========

Runs the entire raw-reads-to-clean-assembly workflow in a single command: trim -> filter -> assemble -> vecscreen -> sourpurge -> rmdup -> polish -> sort -> assess. See
:doc:`../workflow` for the full diagram and per-step input/output description.

Algorithm
=========

``pipeline`` is an orchestrator, not a separate implementation -- it calls each submodule's
``run()`` function directly, in-process. Every step starts from exactly the defaults it has on its
own command line (read from the ``AAFTF <step>`` parser, so ``pipeline`` cannot drift from them);
only the options given to ``pipeline`` and the file names that chain the steps together override
them, plus ``pipe=True`` so intermediate steps don't print "next command" hints.

Before each step, ``pipeline`` checks whether that step's expected output file already exists; if
so, the step is **skipped** (with a log message) rather than re-run. This makes an interrupted
pipeline run resumable: simply re-invoke the same ``AAFTF pipeline`` command and it will pick up
after the last completed step. After each step, output existence is re-checked and the pipeline
aborts with an error if the expected output was not produced.

Steps intentionally **not** included in ``pipeline`` (marked "(Optional)" in ``AAFTF -h``) -- run
these manually if needed: :doc:`mito` (to screen mitochondrial reads out, run it before
:doc:`filter` and pass its output to ``filter --screen_local``), :doc:`fcs_screen`, :doc:`fcs_gx_purge` (alternatives/complements to :doc:`vecscreen` /
:doc:`sourpurge`), :doc:`depth` (coverage QC of the final assembly), :doc:`fix_tbl` (post-FCS
annotation coordinate fixups), and :doc:`database` (run once, ahead of time, to populate
``$AAFTF_DB``).

Cutoffs / defaults
===================

``pipeline`` uses each step's own defaults (:doc:`trim`, :doc:`filter`,
:doc:`assemble`, :doc:`vecscreen`, :doc:`sourpurge`, :doc:`rmdup`, :doc:`polish`, :doc:`sort`,
:doc:`assess`) except where noted:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-m/--memory``
     - each step's own default
     - Passed to every step that has ``-m/--memory`` (trim, filter, assemble, polish)
   * - ``-c/--cpus``, ``-w/--workdir``, ``-v/--verbose``
     - each step's own default
     - Passed to every step that has the option
   * - ``--method``
     - spades
     - Assembler method (spades / megahit / unicycler)
   * - ``-ml/--minlen``
     - 75
     - Minimum read length kept by trim
   * - ``-mc/--mincontiglen``
     - 500
     - Minimum contig length kept by rmdup and by the final sort step
   * - ``--mincovpct``
     - 5
     - sourpurge low-coverage removal threshold (percent of N50-contig coverage)

Invocation
==========

.. code-block:: text

    AAFTF pipeline -1 READ1 [-2 READ2] -o BASENAME -p PHYLUM [PHYLUM ...]
                   [-c CPUS] [-m MEMORY] [-ml MINLEN]
                   [-mc MINCONTIGLEN] [--method {spades,megahit,unicycler}]
                   [-a ACCESSIONS ...] [-u URLS ...] [--sourdb PATH]
                   [--mincovpct PCT] [-w WORKDIR]
                   [--assembler_args ARG ...] [--tmpdir DIR] [-v] [--pipe]

``-1/--read1``, ``-o/--out`` (basename), and ``-p/--phylum`` are required.

Example
=======

.. code-block:: bash

    export AAFTF_DB=~/lib/AAFTF_DB

    AAFTF pipeline \
        -1 reads/STRAINX_R1.fq.gz -2 reads/STRAINX_R2.fq.gz \
        -o STRAINX -c 24 -m 96 \
        --phylum Ascomycota

This produces, in sequence: ``STRAINX_1P.fastq.gz``/``STRAINX_2P.fastq.gz`` (trim),
``STRAINX_filtered_1.fastq.gz``/``_2.fastq.gz``
(filter), ``STRAINX.spades.fasta`` (assemble), ``STRAINX.vecscreen.fasta`` (vecscreen),
``STRAINX.sourpurge.fasta`` (sourpurge), ``STRAINX.rmdup.fasta`` (rmdup),
``STRAINX.polish.fasta`` (polish), and finally ``STRAINX.final.fasta`` with printed/``assess``
summary statistics.

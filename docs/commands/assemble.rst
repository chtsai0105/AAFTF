========
assemble
========

Runs a de novo genome assembler on cleaned (trimmed + filtered) reads.

Algorithm
=========

The assembler is selected with ``--method``:

* **spades** (default) -- ``spades.py --mem <GB>`` in ``--isolate`` mode (default; disable with
  ``--no-isolate``). SPAdes does not allow ``--careful`` together with ``--isolate``, so
  ``--careful`` is only added with ``--no-isolate`` (and dropped with ``--no-careful``).
  ``--cov-cutoff auto`` is added unless ``--meta`` is passed via ``--assembler_args``; merged/single
  reads from ``--merged`` are included. If the
  ``--workdir`` already exists from a prior run, AAFTF instead resumes with
  ``spades.py --restart-from last`` rather than restarting from scratch. Final assembly is copied
  from SPAdes' ``scaffolds.fasta``.
* **megahit** -- fast De Bruijn graph assembler, generally lower accuracy than SPAdes but much
  faster/lower memory; ``-m`` GB is converted to bytes for ``megahit --memory``. Does not support
  resuming: if the output folder already exists, AAFTF stops with an error before running MEGAHIT
  (remove it or pass another ``-w``).
  Final assembly is copied from ``final.contigs.fa``.
* **unicycler** -- wraps SPAdes with additional scaffolding logic; supports combining
  short reads with ``--longreads`` (hybrid assembly). With paired reads, ``--merged`` reads are
  passed as ``--unpaired`` alongside the pairs (ignored for single-end input). Final assembly is
  copied from ``assembly.fasta``.

If the assembler's expected output file is missing after the run, a ``RuntimeError`` is raised
(exit code 1) pointing at the assembler log in the work directory.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``--method``
     - spades
     - spades / megahit / unicycler
   * - ``-m/--memory``
     - 32 (GB)
     - Integer GB; SPAdes ``--mem`` / MEGAHIT ``--memory`` (converted to bytes)
   * - ``--no-isolate``
     - isolate on
     - SPAdes ``--isolate`` mode (recommended for high-coverage isolate data)
   * - ``--no-careful``
     - careful on
     - SPAdes ``--careful``; only applied when ``--no-isolate`` is given
   * - SPAdes coverage cutoff
     - ``auto``
     - ``--cov-cutoff auto`` unless running in ``--meta`` mode

Invocation
==========

.. code-block:: text

    AAFTF assemble -1 FASTQ -o FASTA [-2 FASTQ] [--method {spades,megahit,unicycler}]
                   [-w DIR] [-c INT] [-m GB] [--merged FASTQ] [-lr FASTQ]
                   [--no-careful] [--no-isolate]
                   [--tmpdir DIR] [--assembler_args ARG]
                   [-q] [-v]

``-1/--read1`` and ``-o/--out`` (output assembly FASTA) are required. ``--assembler_args`` takes
one argument and may be repeated to pass several raw arguments to the assembler. ``-lr/--longreads``
is used by Unicycler only. The step log is written to ``<workdir>/assemble.log``.

Example
=======

.. code-block:: bash

    AAFTF assemble -c 24 --memory 96 \
        --read1 reads_filtered_1.fastq.gz --read2 reads_filtered_2.fastq.gz \
        -o genomes/STRAINX.spades.fasta -w working_AAFTF/spades_STRAINX

    # Hybrid short+long read assembly with Unicycler
    AAFTF assemble --method unicycler -c 24 \
        --read1 reads_filtered_1.fastq.gz --read2 reads_filtered_2.fastq.gz \
        --longreads ont_reads.fastq.gz \
        -o genomes/STRAINX.unicycler.fasta

Next step: :doc:`vecscreen`.

=========
sourpurge
=========

Purges contaminant contigs from an assembly using two independent signals: sourmash-based
taxonomic classification, and (if reads are provided) unusually low read coverage relative to the
bulk of the assembly.

Algorithm
=========

1. **Taxonomy screen.** Computes a per-contig sourmash sketch (``sourmash compute -k KMER
   --scaled=1000 --singleton``) and classifies each contig with
   ``sourmash lca classify`` against a reference LCA database. The database is either
   user-supplied (``--sourdb``) or the one installed by ``AAFTF database`` for ``--sourdb_type``:
   ``gbk`` -> ``sm_gbk`` (GenBank k=31, the default), ``gtdbrep`` -> ``sm_gtdbrep``, ``gtdb`` ->
   ``sm_gtdb`` (exit code 2 if it is not installed). Only rows with status exactly ``found`` are
   classified: a contig is marked for removal if none of the names given to ``-p/--phylum``
   matches (exactly, case-sensitive) any rank of its lineage. Contigs with *no* classification
   (``nomatch``) are kept (sourmash can't classify novel/divergent sequence).
2. **Low-coverage screen** (only run if ``-1/--read1`` reads are supplied). Maps reads to the
   taxonomy-filtered assembly with ``bwa mem`` + ``samtools sort``, computes per-contig read
   coverage with ``samtools bedcov``, then computes the mean coverage of contigs at/above the
   assembly's **N50** contig length as a "true genome coverage" baseline. Any contig with
   coverage <= ``mincovpct%`` of that N50-contig average coverage is marked for removal (default
   5% -- i.e. contigs covered at less than 1/20th of the bulk-genome sequencing depth, typically
   evidence of a low-level contaminant or index-hopping/cross-contamination artifact rather than
   the target organism). If no N50-length contig has coverage data, the coverage filter is
   skipped with a warning.
3. Contigs flagged by either screen are dropped; everything else is written to ``-o/--outfile``.
   A copy of the raw sourmash classification CSV is saved as
   ``{input_basename}.sourmash-taxonomy.csv`` next to the input file for review.

Intermediate files and ``sourpurge.log`` are written to the work directory, which is kept when
given with ``-w``, with ``-v``, or when the run fails. Unless ``-q`` is given, the suggested next
command (:doc:`rmdup`) is printed.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``-p/--phylum``
     - *(required)*
     - Phylum/phyla to **keep** (e.g. ``Ascomycota``); anything classified outside this set is
       purged
   * - ``-k/--kmer``
     - 31
     - sourmash k-mer size (must match the database's k-mer size)
   * - ``--sourdb_type``
     - gbk
     - gbk (GenBank) / gtdbrep / gtdb sourmash LCA database from ``AAFTF database``
       (ignored when ``--sourdb`` is given)
   * - ``-mc/--mincovpct``
     - 5
     - Minimum coverage, as a percent of the N50-contig average coverage, to avoid low-coverage
       removal

Invocation
==========

.. code-block:: text

    AAFTF sourpurge -i FASTA -o FASTA -p PHYLUM [PHYLUM ...] [-w DIR]
                    [-1 FASTQ] [-2 FASTQ] [--sourdb FILE] [-k KMER] [-mc PCT]
                    [--sourdb_type {gbk,gtdbrep,gtdb}] [--just-show-taxonomy]
                    [-c INT] [-q] [-v]

``-i/--input``, ``-o/--outfile``, and ``-p/--phylum`` are required. Providing ``-1/--read1`` (and
``-2/--read2`` for paired data) enables the low-coverage screen; without reads, only the taxonomy
screen runs. ``--just-show-taxonomy`` logs the distinct lineages found, saves
``{input_basename}.sourmash-taxonomy.csv`` next to the input, and returns without purging or writing
an output assembly (useful to sanity-check what phyla are present before choosing ``--phylum``).

Example
=======

.. code-block:: bash

    AAFTF sourpurge -c 24 -i genomes/STRAINX.vecscreen.fasta -o genomes/STRAINX.sourpurge.fasta \
        --read1 reads_filtered_1.fastq.gz --read2 reads_filtered_2.fastq.gz \
        --phylum Ascomycota

    # Preview classifications before deciding which phylum to keep
    AAFTF sourpurge -i genomes/STRAINX.vecscreen.fasta -o /dev/null \
        --phylum Ascomycota --just-show-taxonomy

Next step: :doc:`rmdup`.

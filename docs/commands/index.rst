=================
Command Reference
=================

.. toctree::
   :maxdepth: 1

   download
   trim
   mito
   filter
   assemble
   vecscreen
   fcs_screen
   fcs_gx_purge
   sourpurge
   rmdup
   polish
   sort
   assess
   depth
   fix_tbl
   pipeline

.. list-table::
   :header-rows: 1
   :widths: 20 16 46

   * - Canonical name
     - Key tools
     - Purpose
   * - :doc:`download`
     - (urllib only)
     - Fetch/cache reference databases into ``$AAFTF_DB``
   * - :doc:`trim`
     - bbduk.sh, trimmomatic, fastp
     - Adapter + quality trimming
   * - :doc:`mito`
     - NOVOPlasty, minimap2
     - De novo mitochondrial genome assembly
   * - :doc:`filter`
     - bbduk.sh, bowtie2, bwa, minimap2, samtools
     - Remove contaminant/PhiX reads
   * - :doc:`assemble`
     - spades.py, megahit, unicycler
     - Assemble cleaned reads
   * - :doc:`vecscreen`
     - blastn, makeblastdb
     - Vector/Euk/Prok/Mito contamination screen
   * - :doc:`fcs_screen`
     - run_fcsadaptor.sh (singularity/docker)
     - NCBI FCS-adaptor vector trimming
   * - :doc:`fcs_gx_purge`
     - run_gx.py (NCBI FCS-GX)
     - NCBI FCS-GX genomic cross-contamination purge
   * - :doc:`sourpurge`
     - sourmash, bwa, samtools
     - Taxonomy + low-coverage contig purge
   * - :doc:`rmdup`
     - minimap2
     - Remove redundant/duplicate contigs
   * - :doc:`polish`
     - polypolish, pypolca, nextPolish2, racon, bwa, samtools, minimap2, yak, freebayes
     - Short/long-read assembly polishing
   * - :doc:`sort`
     - (BioPython only)
     - Sort contigs by length, rename headers
   * - :doc:`assess`
     - (BioPython only)
     - Assembly completeness / summary statistics
   * - :doc:`depth`
     - minimap2, bwa, samtools, mosdepth
     - Per-contig read-depth report + outlier flags
   * - :doc:`fix_tbl`
     - (none)
     - Adjust NCBI ``.tbl`` feature coordinates after FCS trimming
   * - :doc:`pipeline`
     - all of the above (except fcs_screen/fcs_gx_purge/depth)
     - Run the full trim -> ... -> assess pipeline in one command

Every subcommand accepts ``-v/--debug`` (verbose logging; also usually retains temp working
directories) and ``--pipe`` (suppress the "your next command might be" hint printed at the end --
set automatically when a step is invoked from inside ``AAFTF pipeline``).

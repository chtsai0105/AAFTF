====
mito
====

De novo assembles the mitochondrial genome from trimmed paired-end Illumina reads, independently
of the nuclear genome assembly, using NOVOPlasty's seed-and-extend organelle assembler.

Algorithm
=========

1. Estimates read length from the input FASTQ (``estimate_read_length``).
2. Selects the **seed**, the sequence NOVOPlasty starts assembling from and extends with matching
   reads: a user-supplied ``--seed`` FASTA (e.g. a gene, or the mitochondrial genome of a
   related species) or, by default, a bundled *Aspergillus nidulans* cytochrome-b (cob) fragment
   (``aaftf/data/mito-seed.fasta``).
3. Writes a NOVOPlasty config file from the ``aaftf/data/novoplasty-config.txt`` template
   (substituting project name, min/max genome length, max memory, seed, read length, and the
   forward/reverse FASTQ paths) and runs ``NOVOPlasty.pl -c novo-config.txt``.
4. Parses NOVOPlasty's output directory for (in priority order) a
   ``Circularized_assembly_*``, ``Contigs_1_*``, or ``Uncircularized_assemblies_*`` file.
5. If circularized, rotates/reorients the genome to start at a chosen gene: aligns a start
   sequence (``--starting``; by default a bundled spoa consensus of fungal cytochrome-b genes,
   ``aaftf/data/mito-start-cob.fasta``) against the assembly with
   ``minimap2 -x map-ont``, and rotates the sequence to begin at that alignment's start
   coordinate (reverse-complementing if the hit is on the minus strand). Rotation is skipped
   (contig kept as-is, with a warning) if zero or multiple alignments are found, or if the
   computed rotation offset falls outside the sequence length.
6. If not circularized, all contigs are renumbered and concatenated into the output FASTA as-is,
   with a warning that circularization failed.

Requires ``NOVOPlasty.pl`` and ``minimap2`` on ``$PATH``; exits immediately if either is missing.

Cutoffs / defaults
===================

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Meaning
   * - ``--minlen``
     - 10000
     - Minimum expected mitochondrial genome size (NOVOPlasty search bound)
   * - ``--maxlen``
     - 100000
     - Maximum expected mitochondrial genome size (NOVOPlasty search bound)
   * - ``-m/--memory``
     - 8 (GB)
     - NOVOPlasty max RAM, passed into the generated config
   * - ``--seed``
     - bundled *A. nidulans* cob fragment (``aaftf/data/mito-seed.fasta``)
     - Where NOVOPlasty starts assembling (the sequence it extends)
   * - ``--starting``
     - bundled consensus of fungal cob genes (``aaftf/data/mito-start-cob.fasta``)
     - Start gene the finished circular genome is rotated to begin at (no effect on assembly)

Invocation
==========

.. code-block:: text

    AAFTF mito -l LEFT -r RIGHT [-o OUT] [--minlen N] [--maxlen N]
               [-s/--seed FASTA] [--starting FASTA] [-m MEMORY]
               [-w WORKDIR] [-q] [-v] [--pipe]

``-l/--left`` and ``-r/--right`` are required (``mito`` only supports paired-end data); ``-o/--out``
defaults to ``mito.fasta``.

Example
=======

.. code-block:: bash

    AAFTF mito -l reads_trimmed/STRAINX_1P.fastq.gz -r reads_trimmed/STRAINX_2P.fastq.gz \
        -o STRAINX.mito.fasta

    # Seed with a related species' mitogenome, and rotate the result to start at its nad1 gene
    AAFTF mito -l STRAINX_1P.fastq.gz -r STRAINX_2P.fastq.gz -o STRAINX.mito.fasta \
        --seed related_species_mito.fasta --starting related_species_nad1.fasta

``mito`` is not part of ``pipeline``. To keep mitochondrial reads out of the nuclear assembly,
run it before :doc:`filter` and pass its output as ``filter --screen_local STRAINX.mito.fasta``.

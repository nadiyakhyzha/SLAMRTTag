# SLAMRTTag
Part 1 – RT&Tag Pipeline: Defining Transcripts Enriched in a Nuclear Compartment

This pipeline processes unlabeled RT&Tag data to identify transcripts enriched near a nuclear compartment of interest. Reads are aligned to the genome, quantified at the transcript level, and subjected to differential enrichment analysis.

    Inputs:

        FASTQ files

        Genome assembly (FASTA)

        Gene annotation (GTF)

    Outputs:

        Aligned BAM files

        BigWig coverage tracks

        FeatureCounts tabular read counts

        List of differentially enriched transcripts

Part 2 – SLAM-RT&Tag Pipeline: Measuring Transcript Dynamics

This pipeline analyzes s4U-labeled and unlabeled SLAM-RT&Tag samples to assess transcript dynamics within nuclear compartments. Reads are aligned using the SLAM-DUNK pipeline, which detects T>C conversions. Resulting BAM files are processed with bam2bakR to generate browser tracks (tdr files) with T>C conversion mark-up, followed by BakR analysis to perform differential kinetics analysis.

    Inputs:

        FASTQ files

        Genome assembly (FASTA)

        Gene annotation (GTF)

    Outputs:

        Aligned BAM files

        Summary of T>C conversion rates

        TDR track files

        cB file for BakR

        Differential kinetic analysis results

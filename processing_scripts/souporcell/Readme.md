# Souporcell README

1. `demux_prep_vcfs.sh` - Merge in filter variants by batch
2. `mpileup.sh` - quantify RNA reads per variant sites
     - input: .fa file, .txt file of .cram file path, .vcf file (zipped)
     - output: mpileup file
2. `mpileup_plots.R` - visualise the mpileup output and make .bed files
     - input: mpileup file
     - output: filtered .bed file
3. `mpileup_vcf_filter.sh` - filter the VCFs using the .bed files
     - input: filtered .bed file, .vcf file to be filtered (zipped)
     - output: filtered .vcf file (zipped)
4. `cramtobam.sh` - convert .cram file to .bam for souporcell
     - input:.fa file path, .cram file
     - output: .bam file
5. `run_souporcell.sh` - running souporcell for demultiplexing
     - input: .bam file, .fa file, output file path, no. of samples, filtered .vcf file (unzipped), sample names
     - output: ambient_rna.txt - predicted percentage of ambient RNA; clusters.tsv - file with droplet type and cluster assignment; Genotype_ID_key.txt - cluster correlation to sample; ref_clust_pearson_correlation.png - heatmap of cluster-sample correlation; ref_clust_pearson_correlations.tsv - table of cluster-sample correlation; results.tsv - summary file with droplet type and cluster-sample assignment

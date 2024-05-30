#!/bin/bash

# keeping only snRNAs 
sed '$d' data/ensembl/brainvar_rRNA_snRNA_tpm5.txt | sed 's/ /\t/g' > data/ensembl/brainvar_snRNA_tpm5.txt

# making bed format
awk '{print $1 "\t" $2 "\t" $3}' data/ensembl/brainvar_snRNA_tpm5.txt | tail -n+2 > data/ensembl/brainvar_snRNA_tpm5.bed

# retrieving variants within snRNA regions
bcftools view data/ukb/UKB_vars_500k.vcf.gz -R data/ensembl/brainvar_snRNA_tpm5.bed -Oz -o data/ukb/UKB_vars_snRNA.vcf.gz

# calculating normalised proportion of observed SNVs for snRNA variants
Rscript calculate_depletion.R 'data/ukb/UKB_vars_snRNA.vcf.gz' 'data/ensembl/brainvar_snRNA_tpm5.txt' 'chr' 'start' 'end' 'gene_id' 18 'data/prop_observed_snvs/snrnas.txt'

# making bed file for intergenic regions
tail -n+2 data/ensembl/hg38_ncbiRefSeq_chr12_intergenic_nonoverlap_randseq_no_centremere_updated.txt| awk  -F' ' -v OFS='\t' '{print $1, $4, $5}' > data/ensembl/intergenic_regions.bed

# retrieving variants within intergenic regions
bcftools view data/ukb/UKB_vars_500k.vcf.gz -R data/ensembl/intergenic_regions.bed -Oz -o data/ukb/UKB_vars_chr12_intergenic.vcf.gz

# add region number to intergenic regions file for processing
awk -v OFS='\t' 'BEGIN {print "chrom", "start", "end", "region_id"} {print $0, "region" NR}' data/ensembl/intergenic_regions.bed > 'data/ensembl/intergenic_regions.txt'

# calculating normalised proportion of observed SNVs for intergenic region variants
Rscript calculate_depletion.R 'data/ukb/UKB_vars_chr12_intergenic.vcf.gz' 'data/ensembl/intergenic_regions.txt' 'chrom' 'start' 'end' 'region_id' 18 'data/prop_observed_snvs/intergenic_chr12.txt'

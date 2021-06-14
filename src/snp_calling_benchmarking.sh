#SHOULD BE RAN AFTER RUNNING THE src/methylation_calling_benchmarking.sh script. As this contains the alignment of the denovo 
#cluster to the reference genome!!! If not done in this order the denovo and legacy branch benchmarking won't work.

#Create a conda environment which contains all the dependencies for running the SNP calling benchmarking
conda create env -f src/env/benchmarking.yaml

#SNP Calling
#Download ascencions vcf files
#Col-0
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf_with_quality_reference/6909_snp_short_indel_with_quality_reference.vcf.gz
mv 6909_snp_short_indel_with_quality_reference.vcf.gz data/snp-calls/Col.vcf.gz
#Duplicate columns with sample names to make it easy to use vcfeval to asses SNP calling performance
paste <(zcat data/snp-calls/Col.vcf.gz) \
<(zcat data/snp-calls/Col.vcf.gz | cut -f10 | sed 's/6909/Col0_10/') \
<(zcat data/snp-calls/Col.vcf.gz | cut -f10 | sed 's/6909/Col0_14/') \
<(zcat data/snp-calls/Col.vcf.gz | cut -f10 | sed 's/6909/Col0_1/') \
<(zcat data/snp-calls/Col.vcf.gz | cut -f10 | sed 's/6909/Col0_23/') \
<(zcat data/snp-calls/Col.vcf.gz | cut -f10 | sed 's/6909/Col0_27/') \
<(zcat data/snp-calls/Col.vcf.gz | cut -f10 | sed 's/6909/Col0_36/') \
<(zcat data/snp-calls/Col.vcf.gz | cut -f10 | sed 's/6909/Col0_44/') \
 <(zcat data/snp-calls/Col.vcf.gz | cut -f10 | sed 's/6909/Col0_53/') | sed "s|\t##.*$||" | bgzip -c > data/snp-calls/ColRef.vcf.gz
bcftools tabix data/snp-calls/ColRef.vcf.gz

#Ler0
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf_with_quality_reference/7213_snp_short_indel_with_quality_reference.vcf.gz
mv 7213_snp_short_indel_with_quality_reference.vcf.gz data/snp-calls/Ler.vcf.gz
paste <(zcat data/snp-calls/Ler.vcf.gz) \
<(zcat data/snp-calls/Ler.vcf.gz | cut -f10 | sed 's/7213/Ler0_21/') \
<(zcat data/snp-calls/Ler.vcf.gz | cut -f10 | sed 's/7213/Ler0_32/') \
<(zcat data/snp-calls/Ler.vcf.gz | cut -f10 | sed 's/7213/Ler0_5/') | sed "s|\t##.*$||" | bgzip -c > data/snp-calls/LerRef.vcf.gz
bcftools tabix data/snp-calls/LerRef.vcf.gz

#Cvi 
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf_with_quality_reference/6911_snp_short_indel_with_quality_reference.vcf.gz
mv 6911_snp_short_indel_with_quality_reference.vcf.gz data/snp-calls/Cvi.vcf.gz
paste <(zcat data/snp-calls/Cvi.vcf.gz) \
<(zcat data/snp-calls/Cvi.vcf.gz | cut -f10 | sed 's/6911/Cvi0_11/') \
<(zcat data/snp-calls/Cvi.vcf.gz | cut -f10 | sed 's/6911/Cvi0_16/') \
<(zcat data/snp-calls/Cvi.vcf.gz | cut -f10 | sed 's/6911/Cvi0_20/') \
<(zcat data/snp-calls/Cvi.vcf.gz | cut -f10 | sed 's/6911/Cvi0_31/') \
<(zcat data/snp-calls/Cvi.vcf.gz | cut -f10 | sed 's/6911/Cvi0_39/') \
<(zcat data/snp-calls/Cvi.vcf.gz | cut -f10 | sed 's/6911/Cvi0_3/') \
<(zcat data/snp-calls/Cvi.vcf.gz | cut -f10 | sed 's/6911/Cvi0_48/') \
<(zcat data/snp-calls/Cvi.vcf.gz | cut -f10 | sed 's/6911/Cvi0_49/') | sed "s|\t##.*$||" | bgzip -c > data/snp-calls/CviRef.vcf.gz
bcftools tabix data/snp-calls/CviRef.vcf.gz

#Ei2
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf_with_quality_reference/6915_snp_short_indel_with_quality_reference.vcf.gz
mv 6915_snp_short_indel_with_quality_reference.vcf.gz data/snp-calls/Ei2.vcf.gz
paste <(zcat data/snp-calls/Ei2.vcf.gz) \
<(zcat data/snp-calls/Ei2.vcf.gz | cut -f10 | sed 's/6915/Ei2_18/') \
<(zcat data/snp-calls/Ei2.vcf.gz | cut -f10 | sed 's/6915/Ei2_19/') \
<(zcat data/snp-calls/Ei2.vcf.gz | cut -f10 | sed 's/6915/Ei2_2/') \
<(zcat data/snp-calls/Ei2.vcf.gz | cut -f10 | sed 's/6915/Ei2_34/') \
<(zcat data/snp-calls/Ei2.vcf.gz | cut -f10 | sed 's/6915/Ei2_38/') \
<(zcat data/snp-calls/Ei2.vcf.gz | cut -f10 | sed 's/6915/Ei2_51/') \
<(zcat data/snp-calls/Ei2.vcf.gz | cut -f10 | sed 's/6915/Ei2_9/') | sed "s|\t##.*$||" | bgzip -c > data/snp-calls/Ei2Ref.vcf.gz
bcftools tabix data/snp-calls/Ei2Ref.vcf.gz

#Gu0
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf_with_quality_reference/6922_snp_short_indel_with_quality_reference.vcf.gz
mv 6922_snp_short_indel_with_quality_reference.vcf.gz data/snp-calls/Gu0.vcf.gz
paste <(zcat data/snp-calls/Gu0.vcf.gz) \
<(zcat data/snp-calls/Gu0.vcf.gz | cut -f10 | sed 's/6922/Gu0_15/') \
<(zcat data/snp-calls/Gu0.vcf.gz | cut -f10 | sed 's/6922/Gu0_22/') \
<(zcat data/snp-calls/Gu0.vcf.gz | cut -f10 | sed 's/6922/Gu0_25/') \
<(zcat data/snp-calls/Gu0.vcf.gz | cut -f10 | sed 's/6922/Gu0_35/') \
<(zcat data/snp-calls/Gu0.vcf.gz | cut -f10 | sed 's/6922/Gu0_43/') \
<(zcat data/snp-calls/Gu0.vcf.gz | cut -f10 | sed 's/6922/Gu0_52/') \
<(zcat data/snp-calls/Gu0.vcf.gz | cut -f10 | sed 's/6922/Gu0_6/') \
<(zcat data/snp-calls/Gu0.vcf.gz | cut -f10 | sed 's/6922/Gu0_8/') | sed "s|\t##.*$||" | bgzip -c > data/snp-calls/Gu0Ref.vcf.gz
bcftools tabix data/snp-calls/Gu0Ref.vcf.gz

bcftools tabix data/snp-calls/*Ref.vcf.gz
bcftools merge data/snp-calls/*Ref.vcf.gz data/snp-calls/refAll.vcf.gz
#For pooled do:
cat data/snp-calls/refAll.vcf.gz | sed 's/6909/Col0/' | sed 's/7213/Ler0/' | sed 's/6911/Cvi0/' | sed 's/6915/Ei2/' | sed 's/6922/Gu0/' | bgzip -c  > data/snp-calls/RefAllPooledIncl.vcf

#############################
#SNP_calling reference
#Adjust configRef.yaml input_dir parameter to point to the location of the epiGBS-reference branch output and reference location
#Run snakemake which normalises, filters and performs the evaluation of the snp callign using rtg vcfeval 
snakemake -j 15 -s referenceSNP_calling.smk --use-conda

################################### 
#SNP_calling denovo
#Remove the C24 samples for the VCF file (no 1001 genomes data available)
bcftools query -l /mnt/nfs/bioinfdata/home/NIOO/maartenp/adam/adam-denovo/output/snp_calling/snp.vcf.gz | grep -v "C24" > scratch/samplesDenovo.txt
#Create header from reference to add to the lifted over vcf file 
bcftools view -h /mnt/nfs/bioinfdata/home/NIOO/maartenp/adam/adam-ref/output/snp_calling/snp.vcf.gz > scratch/headerRef.vcf

#Lift over vcf file to the reference coordinates
cat scratch/headerRef.vcf <(python src/liftover_vcf.py scratch/consensus.bam /mnt/nfs/bioinfdata/home/NIOO/maartenp/adam/adam-denovo/output/snp_calling/snp.vcf.gz | cut -f1,2,7-) | bcftools sort | bgzip > results/liftOverAdam.vcf.gz
tabix scratch/liftoverPooled.vcf.gz -f

bcftools norm --check-ref s -f /mnt/nfs/bioinfdata/home/NIOO/maartenp/epiGBS-ref/input/TAIR_10_chr.fa results/liftOverAdam.vcf.gz | bgzip -c  > scratch/denovo/normalisedAdam.vcf.gz
snakemake -s denovoSNP_calling.smk --use-conda -j 16


rtg rocplot --png results/snp_calling/figures/CviDenovo.png  results/snp_calling/denovo/RTGGQ/Cvi*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/Gu0Denovo.png  results/snp_calling/denovo/RTGGQ/Gu0*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/Ei2Denovo.png  results/snp_calling/denovo/RTGGQ/Ei2*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/LerDenovo.png  results/snp_calling/denovo/RTGGQ/Ler*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/ColDenovo.png  results/snp_calling/denovo/RTGGQ/Col*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100

rtg rocplot --png results/snp_calling/figures/CvidenovoNoSquash.png  results/snp_calling/denovoNoSquash/RTGGQ/Cvi*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/Gu0denovoNoSquash.png  results/snp_calling/denovoNoSquash/RTGGQ/Gu0*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/Ei2denovoNoSquash.png  results/snp_calling/denovoNoSquash/RTGGQ/Ei2*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/LerdenovoNoSquash.png  results/snp_calling/denovoNoSquash/RTGGQ/Ler*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/ColdenovoNoSquash.png  results/snp_calling/denovoNoSquash/RTGGQ/Col*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100

#####################

#SNP_calling legacy
#Remove the C24 samples for the VCF file (no 1001 genomes data available)
bcftools query -l /mnt/nfs/bioinfdata/home/NIOO/maartenp/epiGBS-legacy/snp_calling/snp.vcf.gz | grep -v "C24" > scratch/samplesDenovo.txt
#Create header from reference to add to the lifted over vcf file 
bcftools view -h /mnt/nfs/bioinfdata/home/NIOO/maartenp/adam/adam-ref/output/snp_calling/snp.vcf.gz > scratch/headerRef.vcf

#Lift over vcf file to the reference coordinates
cat scratch/headerRef.vcf <(python src/liftover_vcf.py scratch/consensus.bam /mnt/nfs/bioinfdata/home/NIOO/maartenp/epiGBS-legacy/snp_calling/snp_calling/snp.vcf.gz | cut -f1,2,7-) | bcftools sort | bgzip > results/liftOverAdam.vcf.gz
tabix scratch/liftoverPooled.vcf.gz -f

bcftools norm --check-ref s -f /mnt/nfs/bioinfdata/home/NIOO/maartenp/epiGBS-ref/input/TAIR_10_chr.fa results/liftOverAdam.vcf.gz | bgzip -c  > scratch/denovo/normalisedAdam.vcf.gz
snakemake -s denovoSNP_calling.smk --use-conda -j 16



rtg rocplot --png results/snp_calling/figures/CvireferenceFilt.png  results/snp_calling/refFilt/Squash*/Cvi*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/Gu0referenceFilt.png  results/snp_calling/refFilt/Squash*/Gu0*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/Ei2referenceFilt.png  results/snp_calling/refFilt/Squash*/Ei2*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/LerreferenceFilt.png  results/snp_calling/refFilt/Squash*/Ler0*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100
rtg rocplot --png results/snp_calling/figures/ColreferenceFilt.png  results/snp_calling/refFilt/Squash*/Col*/weighted_roc.tsv.gz --precision-sensitivity --zoom=0,0,100,100

Rscript SNP_calling_plots.R

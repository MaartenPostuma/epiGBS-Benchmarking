
bcftools mpileup /home/NIOO/maartenp/epiGBS-denovo/tmp/alignment/Ler0_21_sorted.bam \
--fasta-ref /home/NIOO/maartenp/epiGBS-denovo/output/output_denovo/NNNNref/ref.fa -r 1 | \
grep "#" | sed "s/4.2/4.1/" | grep -v "##FORMAT" | grep -v "##INFO" > input/snpdb.vcf

bis-snp -R ../epiGBS-denovo/output/output_denovo/NNNNref/ref.fa \
        -I ../epiGBS-denovo/tmp/alignment/Ler0_21_sorted.bam \
        -T BisulfiteCountCovariates -knownSites input/snpdb.vcf \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate\
        -recalFile output/Ler0-21_recal_before.csv \
        -nt 1

bis-snp -R ../epiGBS-denovo/output/output_denovo/NNNNref/ref.fa   \
        -I ../epiGBS-denovo/tmp/alignment/Ler0_21_sorted.bam \
        -o output/Ler0_21_recal.bam  -knownSites input/snpdb.vcf\
        -T BisulfiteTableRecalibration \
        -recalFile output/Ler0-21_recal_before.csv  \
        -maxQ 40

bis-snp -R ../epiGBS-denovo/output/output_denovo/NNNNref/ref.fa \
        -T BisulfiteGenotyper \
        -I output/Ler0_21_recal.bam \
        -vfn1 output/Ler-denovo-021.vcf \
        -out_modes EMIT_ALL_CONFIDENT_SITES \
        -stand_call_conf 4\
        -toCoverage 1000 \
        -mmq 20 \
        -mbq 20 \
        -nt 1




bis-snp -R ../epiGBS-ref/input/TAIR_10_chr.fa \ 
        -I ../epiGBS-ref/tmp/alignment/Ler0_21_sorted.bam \
        -T BisulfiteCountCovariates -knownSites input/dbsnp.vcf \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate\
        -recalFile output/Ler0-21_recal_before.csv \
        -nt 1
bis-snp -R ../epiGBS-ref/input/TAIR_10_chr.fa  \
        -I ../epiGBS-ref/tmp/alignment/Ler0_21_sorted.bam \
        -o output/Ler0_21_recal.bam  \
        -T BisulfiteTableRecalibration \
        -recalFile output/Ler0-21_recal_before.csv  \
        -maxQ 
bis-snp -R ../epiGBS-ref/input/TAIR_10_chr.fa \
        -T BisulfiteGenotyper \
        -I output/Ler0_21_recal.bam \
        -vfn1 output/Ler-021.vcf \
        -out_modes EMIT_ALL_CONFIDENT_SITES \
        -stand_call_conf 4\
        -toCoverage 1000 \
        -mmq 20 \
        -mbq 20 \
        -nt 1
source activate snp_benchmarking
bgzip -c output/Ler-021.vcf > output/Ler-021.vcf.gz
bcftools tabix output/Ler-021.vcf.gz -f

bcftools query -l ../SNP-test/output/Ler-021.vcf.gz | sed "s/SAMPLE/7213/" > scratch/samples.txt

bcftools reheader -s scratch/samples.txt ../SNP-test/output/Ler-021.vcf.gz > ../SNP-test/output/Ler-021-reheader.vcf.gz
bcftools norm -f data/ref/TAIR_10_chr.fa ../SNP-test/output/Ler-021-reheader.vcf.gz > ../SNP-test/output/Ler-021-norm.vcf.gz
bcftools view -s 7213 -V indels,mnps,ref,bnd,other ../SNP-test/output/Ler-021-norm.vcf.gz| awk '$0~"^#" || ($10~"^1/1" || $10~"^1/0" || $10~"^0/1")' | bgzip -c > ../SNP-test/output/Ler-021-filt.vcf.gz
bcftools tabix ../SNP-test/output/Ler-021-filt.vcf.gz -f

rtg vcfeval -b scratch/Ler0_23_ref.vcf.gz -c ../SNP-test/output/Ler-021-filt.vcf.gz -t data/ref/TAIR_10_chr.sdf -o results/RTG/Ler0_21 \
--sample=7213 \
--squash-ploidy \
--vcf-score-field=DP

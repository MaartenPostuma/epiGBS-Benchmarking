#####
# init sym links
ln -s ../../../epiGBS-ref/input/ ./
ln -s ../../../finalBenchmarking/data/snp-calls/Cvi.vcf.gz

# build raw reference from 1001genomes (already filtered for 6911 i.e. Cvi-0)
bcftools reheader -s <(echo -e "Cvi-0_11\n") Cvi.vcf.gz |
bcftools norm -f input/TAIR_10_chr.fa -m- |
bcftools view -V indels,mnps,ref,bnd,other -Oz > 6911.snps.vcf.gz
tabix 6911.snps.vcf.gz
# high-confidence set (positions intersected later by RTG)
bcftools view -e 'GT=="./." || GT=="./1" || GT=="0/."' 6911.snps.vcf.gz |
bedtools merge -i stdin > high_confidence.bed


###Create liftover file
##Make consensus clusters into "reads" for minimap
awk '$0~"^>" {if(NR!=1){printf "\n%s\t",$0}else{printf "%s\t",$0}} $0!~"^>" {printf "%s",$0} END{print null}' /home/NIOO/maartenp/adam/adam-denovo/output/output_denovo/consensus_cluster.renamed.fa | awk '{print $1 >> "R1.fa"; print $1 >> "R2.fa"; R1=substr($2,0,(length($2)/2)-1); R2=substr($2,length($2)/2); print R1 >> "R1.fa"; print R2 >> "R2.fa"}';
paste <(grep "^>" R2.fa) <(grep -v "^>" R2.fa | tr ACTG TGAC | rev) | tr "\t" "\n" > consensus_cluster_with_Ns_2.fa && rm R2.fa;
mv R1.fa consensus_cluster_with_Ns_1.fa
##Run minimap
minimap2 -ax sr input/TAIR_10_chr.fa consensus_cluster_with_Ns_1.fa consensus_cluster_with_Ns_2.fa > consensus_cluster_Ns.sam
samtools view -Sb consensus_cluster_Ns.sam   > consensus.bam
#Calculate depth of the denovo bam
samtools sort ../output/alignment/Cvi0_11_trimmed_filt_merged.1_bismark_bt2_pe.bam > Cvi0_11.bam
samtools index Cvi0_11.bam
samtools depth Cvi0_11.bam | awk '$3>0{print $1,$2-1,$2,$3}' | sed 's/ /\t/g' > Cvi0_11_all.depth  #outputs space delimited instead of tab?


#Run the lift over script
python ../../../finalBenchmarking/src/liftover_bed.py consensus.bam Cvi0_11_all.depth > AllLiftover.depth


# calculate per-strand depth
echo -e "83\n99\n147\n163\n" |
xargs -n1 -P4 -i sh -c "samtools view -bf '{}' Cvi0_11.bam > '{}'.bam"
samtools merge Cvi0_11.CT.bam 99.bam 147.bam &
samtools merge Cvi0_11.GA.bam 83.bam 163.bam

bedtools intersect -a <(samtools depth -Q10 Cvi0_11.CT.bam | awk -v OFS="\t" '$3>0{print $1,$2-1,$2,$3}') \
-b <(samtools depth -Q10 Cvi0_11.GA.bam | awk -v OFS="\t" '$3>0{print $1,$2-1,$2,$3}') -u > Cvi0_11.depth

python ../../../finalBenchmarking/src/liftover_bed.py consensus.bam Cvi0_11.depth > Filtered.depth

# intersect Reference with depth file
cat <(bcftools view -h 6911.snps.vcf.gz) \
<(bedtools intersect -a 6911.snps.vcf.gz -b Filtered.depth -u) |
bgzip -c > Cvi0_11.Ref.Filtered.vcf.gz
tabix Cvi0_11.Ref.Filtered.vcf.gz

# intersect Reference with depth file
cat <(bcftools view -h 6911.snps.vcf.gz) \
<(bedtools intersect -a 6911.snps.vcf.gz -b AllLiftover.depth -u) |
bgzip -c > Cvi0_11.Ref.vcf.gz
tabix Cvi0_11.Ref.vcf.gz

#####
# normalising and filtering
bcftools view -s Cvi0_11 ../output/snp_calling/snp.vcf.gz | bgzip -c > Cvi0_11.vcf.gz
tabix Cvi0_11.vcf.gz
#Lift over vcf file to the reference coordinates

bcftools view -h -s Cvi0_11 /mnt/nfs/bioinfdata/home/NIOO/maartenp/adam/adam-ref/output/snp_calling/snp.vcf.gz > headerRef.vcf

cat headerRef.vcf <(python ../../../finalBenchmarking/src/liftover_vcf.py consensus.bam Cvi0_11.vcf.gz | cut -f1,2,7-) | bcftools sort | bgzip > liftOver.vcf.gz
tabix liftOver.vcf.gz -f

bcftools reheader -s <(echo -e "Cvi-0_11\n") liftOver.vcf.gz |
bcftools norm --check-ref s -f input/TAIR_10_chr.fa -m- |
bcftools view -V indels,mnps,ref,bnd,other -Oz > Cvi0_11_norm.vcf.gz
tabix Cvi0_11_norm.vcf.gz -f
# intersect Call set with depth file
cat <(bcftools view -h Cvi0_11_norm.vcf.gz) \
<(bedtools intersect -a Cvi0_11_norm.vcf.gz -b Filtered.depth -u) |
bgzip -c > Cvi0_11_filtered.vcf.gz
tabix Cvi0_11_filtered.vcf.gz -f

#####
# run RTG vcfeval

rtg vcfeval -b Cvi0_11.Ref.Filtered.vcf.gz -c Cvi0_11_filtered.vcf.gz -t input/TAIR_10_chr.fa.sdf -o DP_Filt \
  --vcf-score-field=DP \
  --evaluation-regions=high_confidence.bed 


rtg vcfeval -b Cvi0_11.Ref.Filtered.vcf.gz -c Cvi0_11_filtered.vcf.gz -t input/TAIR_10_chr.fa.sdf -o DP_Filt_squash \
  --vcf-score-field=DP \
  --evaluation-regions=high_confidence.bed --squash-ploidy

rtg vcfeval -b Cvi0_11.Ref.vcf.gz -c Cvi0_11_norm.vcf.gz -t input/TAIR_10_chr.fa.sdf -o DP_NoFilt \
  --vcf-score-field=DP \
  --evaluation-regions=high_confidence.bed 


rtg vcfeval -b Cvi0_11.Ref.vcf.gz -c Cvi0_11_norm.vcf.gz -t input/TAIR_10_chr.fa.sdf -o DP_NoFilt_squash \
  --vcf-score-field=DP \
  --evaluation-regions=high_confidence.bed --squash-ploidy


  rtg rocplot --png=ref_DP.png --precision-sensitivity --zoom=0,0,100,100 DP*/weighted_roc.tsv.gz

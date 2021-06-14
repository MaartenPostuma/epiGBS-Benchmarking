#Activate conda environment
conda env create -f src/env/epiGBS-benchmarking.yaml
source activate epiGBS-benchmarking.yaml

#Download col-0 data #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43857 for other ascessions
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1085nnn/GSM1085222/suppl/GSM1085222_mC_calls_Col_0.tsv.gz
mv GSM1085222_mC_calls_Col_0.tsv.gz data/methylation-calls/Col0-calls.tsv.gz

#Download C24 data
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2099nnn/GSM2099232/suppl/GSM2099232_allc_6906.tsv.gz
mv GSM2099232_allc_6906.tsv.gz data/methylation-calls/C24-calls.tsv.gz

#Download EI-2 data
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2099nnn/GSM2099234/suppl/GSM2099234_allc_6915.tsv.gz
mv GSM2099234_allc_6915.tsv.gz data/methylation-calls/Ei2-calls.tsv.gz

#Download cvi0
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1085nnn/GSM1085224/suppl/GSM1085224_mC_calls_Cvi_0.tsv.gz
mv GSM1085224_mC_calls_Cvi_0.tsv.gz data/methylation-calls/Cvi0-calls.tsv.gz

#Download Gu0
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1085nnn/GSM1085247/suppl/GSM1085247_mC_calls_Gu_0.tsv.gz
mv GSM1085247_mC_calls_Gu_0.tsv.gz data/methylation-calls/Gu0-calls.tsv.gz

#Download arabidopsis genome
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
mv TAIR10_chr_all.fas data/ref/
sed 's/mitochondria/6/' data/ref/TAIR10_chr_all.fas |sed 's/chloroplast/7/' > data/ref/TAIR_10_chr.fa #remove mitochondria / chloroplast sequences as they break the liftover scripts

#Reference benchmarking

Rscript src/finalBenchMarkingReference.R


###Prepare lift_over data

# process consensus_cluster.renamed.fa ready for PE alignment

awk '$0~"^>" {if(NR!=1){printf "\n%s\t",$0}else{printf "%s\t",$0}} $0!~"^>" {printf "%s",$0} END{print null}' /home/NIOO/maartenp/adam/adam-denovo/output/output_denovo/NNNNref/ref.fa | awk '{print $1 >> "R1.fa"; print $1 >> "R2.fa"; R1=substr($2,0,(length($2)/2)-1); R2=substr($2,length($2)/2); print R1 >> "R1.fa"; print R2 >> "R2.fa"}';
paste <(grep "^>" R2.fa) <(grep -v "^>" R2.fa | tr ACTG TGAC | rev) | tr "\t" "\n" > scratch/consensus_cluster_with_Ns_2.fa && rm R2.fa;
mv R1.fa scratch/consensus_cluster_with_Ns_1.fa


#Align denovo to reference
minimap2 -ax sr data/ref/TAIR_10_chr.fa scratch/consensus_cluster_with_Ns_1.fa scratch/consensus_cluster_with_Ns_2.fa > scratch/consensus_cluster_Ns.sam
samtools view -Sb consensus_cluster_Ns.sam   > consensus.bam


samtools depth Cvi0_11.bam | awk '$3>0{print $1,$2-1,$2,$3}' | sed 's/ /\t/g' > Cvi0_11_all.depth  #outputs space delimited instead of tab?


python ../../../finalBenchmarking/src/liftover_bed.py consensus.bam Cvi0_11_all.depth > AllLiftover.depth

#Run the Rscripts to generate the plots. Make sure to change the file paths to the corresponding file / directory locations!

Rscript src/finalBenchmarkingDenovo.R
Rscript src/finalBenchmarkingReference.R
Rscript src/methylKitPCA.R
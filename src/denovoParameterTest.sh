#denovo paramter test
##This is how it was done for the paper. Since then I've included an option that performs the non denovo/refernece comparisons. 
##If you write paramTest in the config file instead of denovo or reference in the mode field and run the pipeline as normal it
##It will output a figure with the number of clusters, mapping % and % multimapping. For different combinations of the parameters.
##The parameter combinations are in 


git clone epigbs2
for i in {95,97,99};
    do
    for j in {10,50,100};
        do
        cp -r epigbs2/ $i-$j
        cat epigbs2/config.yaml | sed "s/clID/$i/" | sed "s/minDepth/$j/" > $i-$j/config.yaml
        done
    done
#Run epiGBS for a single sample 
snakemake -j 4 /mnt/nfs/bioinfdata/home/NIOO/maartenp/clusterParameterTest/95-10/output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam --use-conda
snakemake -j 4 /mnt/nfs/bioinfdata/home/NIOO/maartenp/clusterParameterTest/95-50/output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam --use-conda
snakemake -j 4 /mnt/nfs/bioinfdata/home/NIOO/maartenp/clusterParameterTest/95-100/output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam --use-conda

snakemake -j 4 /mnt/nfs/bioinfdata/home/NIOO/maartenp/clusterParameterTest/97-10/output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam --use-conda
snakemake -j 4 /mnt/nfs/bioinfdata/home/NIOO/maartenp/clusterParameterTest/97-50/output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam --use-conda
snakemake -j 4 /mnt/nfs/bioinfdata/home/NIOO/maartenp/clusterParameterTest/97-100/output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam --use-conda

snakemake -j 4 /mnt/nfs/bioinfdata/home/NIOO/maartenp/clusterParameterTest/99-10/output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam --use-conda
snakemake -j 4 /mnt/nfs/bioinfdata/home/NIOO/maartenp/clusterParameterTest/99-50/output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam --use-conda
snakemake -j 4 /mnt/nfs/bioinfdata/home/NIOO/maartenp/clusterParameterTest/99-100/output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam --use-conda


for i in {95,97,99};
    do
    for j in {10,50,100};
        do
        awk '$0~"^>" {if(NR!=1){printf "\n%s\t",$0}else{printf "%s\t",$0}} $0!~"^>" {printf "%s",$0} END{print null}' /home/NIOO/maartenp/clusterParameterTest/$i-$j/output/output_denovo/consensus_cluster.renamed.fa |
        awk '{print $1 >> "R1.fa"; print $1 >> "R2.fa"; R1=substr($2,0,(length($2)/2)-1); R2=substr($2,length($2)/2); print R1 >> "R1.fa"; print R2 >> "R2.fa"}';
        paste <(grep "^>" R2.fa) <(grep -v "^>" R2.fa | tr ACTG TGAC | rev) | tr "\t" "\n" > scratch/100consensus_cluster95_with_Ns_2.fa && rm R2.fa;
        mv R1.fa scratch/100consensus_cluster95_with_Ns_1.fa
        #Align denovo to reference
        minimap2 -ax sr data/ref/TAIR_10_chr.fa scratch/100consensus_cluster95_with_Ns_1.fa scratch/100consensus_cluster95_with_Ns_2.fa > scratch/100consensus_cluster95_Ns.sam --secondary=yes
        samtools flagstat scratch/100consensus_cluster95_Ns.sam > results/denovoParameterTest/final$i-$j.txt
        rm scratch/100*
        done
    done


samtools sort output/alignment/Cvi0_3_trimmed_filt_merged.1_bismark_bt2_pe.bam | samtools depth - |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
#Look in the report files and make a csv file containing all the different charateristics:
#"N clusters","Mapping % denovo on ref","No. Clusters map on Ref","Mapping % reads on denovo","% multimapped","Average depth","CG methylation","CHG methylation","CHH methylation"

#
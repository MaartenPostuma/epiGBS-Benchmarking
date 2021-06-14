#Bis-SNP needs some recalibration steps to properly function and output DP values. For this step it needs known snps. 
#But since we are working with a denovo genome these are not available. 
#Luckily bis-snps also accepts empty vcf files as known snps. But this vcf file needs to be vcf4.1. 
#We use mpileup from on the first sorted bam file to generate a functional vcf header.

configfile: "configPooledSNPCalling.yaml"
import pandas as pd
import os
import random
df = pd.read_csv(os.path.join(config["input_dir"],config["barcodes"]), sep='\t', dtype="object").set_index('Sample')
SAMPLES = df.index
genomePlaceHolder = config["genome"].split(".")[0]

rule all:
    input:
        finalVCF=expand("{path}/snp_calling/snp.vcf.gz",path=config["output_dir"]),
        perIndvcfGZ=expand("{path}/snp_calling/{sample}.vcf.gz",path=config["output_dir"],sample=SAMPLES)

rule sort:
    input:
        expand("{indir}/{{sample}}.bam",indir=config["input_dir"])
    output:
        alignment=expand("{out}/alignment/{{sample}}_sorted.bam",out=config["tmp_dir"]),
        aligmentIndex=expand("{out}/alignment/{{sample}}_sorted.bam.bai",out=config["tmp_dir"])
    conda:
        "src/env/samtools.yaml"
    shell:
        """
        samtools sort {input} > {output.alignment}
        samtools index {output.alignment}
        """

rule vcfHeader:
    input:
        alignment=expand("{out}/alignment/{sample}_sorted.bam",out=config["tmp_dir"],sample=SAMPLES[0]),
        aligmentIndex=expand("{out}/alignment/{sample}_sorted.bam.bai",out=config["tmp_dir"],sample=SAMPLES[0]),
        reference=expand("{path}/{genome}",path=config["ref_dir"],genome=config["genome"]),
        AsReferenceDict=expand("{path}/{genome}.dict",path=config["ref_dir"],genome=genomePlaceHolder),
        referenceFai=expand("{path}/{genome}.fai",path=config["ref_dir"],genome=config["genome"])
    output:
        knownSNPs=temp(expand("{out}/alignment/knownsnps.vcf",out=config["tmp_dir"]))
    conda:
        "src/env/bcftools.yaml"
    shell:
        """
        bcftools mpileup {input.alignment} \
        --fasta-ref {input.reference} -r 1 | \
        grep "#" | sed "s/4.2/4.1/" | grep -v "##FORMAT" | grep -v "##INFO" > {output.knownSNPs}
        """

rule countCovariates:
    input:
        alignment=expand("{out}/alignment/{{sample}}_sorted.bam",out=config["tmp_dir"]),
        aligmentIndex=expand("{out}/alignment/{{sample}}_sorted.bam.bai",out=config["tmp_dir"]),
        reference=expand("{path}/{genome}",path=config["ref_dir"],genome=config["genome"]),
        AsReferenceDict=expand("{path}/{genome}.dict",path=config["ref_dir"],genome=genomePlaceHolder),
        referenceFai=expand("{path}/{genome}.fai",path=config["ref_dir"],genome=config["genome"]),
        knownSNPs=expand("{out}/alignment/knownsnps.vcf",out=config["tmp_dir"])
    output:
        recalc=temp(expand("{out}/alignment/{{sample}}_recalc_before.csv",out=config["tmp_dir"]))
    conda:
        "src/env/env_bis-snp0.82.yaml"
    shell:
        """
        bis-snp -R {input.reference} \
        -I {input.alignment} \
        -T BisulfiteCountCovariates -knownSites {input.knownSNPs} \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate\
        -recalFile {output.recalc} \
        -nt 1 
        """
rule recalc:
    input:
        alignment=expand("{out}/alignment/{{sample}}_sorted.bam",out=config["tmp_dir"]),
        aligmentIndex=expand("{out}/alignment/{{sample}}_sorted.bam.bai",out=config["tmp_dir"]),
        reference=expand("{path}/{genome}",path=config["ref_dir"],genome=config["genome"]),
        AsReferenceDict=expand("{path}/{genome}.dict",path=config["ref_dir"],genome=genomePlaceHolder),
        referenceFai=expand("{path}/{genome}.fai",path=config["ref_dir"],genome=config["genome"]),
        knownSNPs=expand("{out}/alignment/knownsnps.vcf",out=config["tmp_dir"]),
        recalc=expand("{out}/alignment/{{sample}}_recalc_before.csv",out=config["tmp_dir"])
    output:
        recalcedBam=temp(expand("{out}/alignment/{{sample}}_recalc.bam",out=config["tmp_dir"]))
    conda:
        "src/env/env_bis-snp0.82.yaml"
    shell:
        """
        bis-snp -R {input.reference}   \
        -I {input.alignment} \
        -o {output.recalcedBam} \
        -T BisulfiteTableRecalibration \
        -recalFile {input.recalc} \
        -maxQ 40 
        """
rule indexRecalc:
    input:
        recalcedBam=expand("{out}/alignment/{{sample}}_recalc.bam",out=config["tmp_dir"])
    output:
        recalcedBai=temp(expand("{out}/alignment/{{sample}}_recalc.bam.bai",out=config["tmp_dir"]))
    conda:
        "src/env/samtools.yaml"
    shell:
        "samtools index {input.recalcedBam}"

# call SNPs from recalibrated alignments
# -XX:ActiveProcessorCount=1  <- ensures not all cores are being used (although still uses more cores then necessarysrc.)
# -out_modes EMIT_ALL_CONFIDENT_SITES (Work in progress (default only outputs SNP sites so when combining you do not find coverage on homozygous ref alleles))

# -stand_call_conf defines the likelihood ratio criteria between best and second best genotype for call to be considered confident
# #TODO Default value is 20 for high depth of coverage. For multiple samples with low coverage (more than 100 samples with 4X coverage), the
# threshold could be defined lower than 10, or even 4. For ultra-high coverage sequencing, such as
# 50X, you could specify higher threshold to obtain higher accuracy. (So maybe play with this value?)

stand_call_conf=config["param_SNPcalling"]["stand_call_conf"]


def getStand_call_conf(stand_call_conf):
    if stand_call_conf == "default" or stand_call_conf == "":
        print(param_mq)
        mq = 20
    else:
        mq = stand_call_conf
    return mq

# -toCoverage 1000 maximum read coverage 
# -toCoverage 1000 maximum read coverage allowed
#-mmq 20 minimal mapping quality
#-mbq 20 minimal base quality
#-nt 1 "Number of threads" actually number of gatk instances started.
rule call_SNPs_reference:
    input:
        alignment=expand("{out}/alignment/{{sample}}_recalc.bam",out=config["tmp_dir"]),
        recalcedBai=expand("{out}/alignment/{{sample}}_recalc.bam.bai",out=config["tmp_dir"]),
        reference=expand("{path}/{genome}",path=config["ref_dir"],genome=config["genome"]),
        AsReferenceDict=expand("{path}/{genome}.dict",path=config["ref_dir"],genome=genomePlaceHolder),
        referenceFai=expand("{path}/{genome}.fai",path=config["ref_dir"],genome=config["genome"])
    output:
        perIndVCF=temp(expand("{path}/snp_calling/{{sample}}_snp.raw.vcf.out",path=config["tmp_dir"])),
        perIndMethStats=temp(expand("{path}/snp_calling/{{sample}}_snp.raw.vcf.out.MethySummarizeList.txt",path=config["tmp_dir"]))
    params:
        stand_call_conf=getStand_call_conf(stand_call_conf)
    log:
        expand("{path}/log/snp_call/{{sample}}_call.log",path=config["output_dir"])
    conda:
        "src/env/bisSNP.yaml"
    threads: 6
    shell:
        "(bis-snp -Xmx4096m -XX:ActiveProcessorCount=1 -T BisulfiteGenotyper -R {input.reference} -I {input.alignment} -vfn1 {output.perIndVCF} -out_modes EMIT_ALL_CONFIDENT_SITES -stand_call_conf {params.stand_call_conf} -toCoverage 1000 -mmq 20 -mbq 20 -nt 1) 2> {log}"

#Combines all different output files
#Order / zip vcf files
#bcftools annotate -x FORMAT/GP Removes a field which causes issues when merging
rule zip_vcfs_reference:
    input:
        perIndVCF=expand("{path}/snp_calling/{{sample}}_snp.raw.vcf.out",path=config["tmp_dir"])
    output:
        perIndvcfGZ=expand("{path}/snp_calling/{{sample}}.vcf.gz",path=config["output_dir"])
    params:
        sample="{sample}"
    conda:
        "src/env/bcftools.yaml"
    shell:
        "bcftools view {input.perIndVCF}| sed 's/SAMPLE/{params.sample}/' | bcftools sort | bcftools annotate -x FORMAT/GP | bgzip -c > {output.perIndvcfGZ}"
#index them so we can use bcftools merge
rule index_vcfs_reference:
    input:
        perIndvcfGZ=expand("{path}/snp_calling/{{sample}}.vcf.gz",path=config["output_dir"])
    output:
        perIndVCFTBI=expand("{path}/snp_calling/{{sample}}.vcf.gz.csi",path=config["output_dir"])
    conda:
        "src/env/bcftools.yaml"
    shell:
        "bcftools index {input.perIndvcfGZ}"

#Merges the individual vcf files
# bcftools view -i 'GT[*]="alt"&&REF!="N"'. Filters site which do not have an alt allele in any individual + sites for which the ref = N
rule merge_vcfs_reference:
    input:
        perIndVCFTBI=expand("{path}/snp_calling/{sample}.vcf.gz.csi",path=config["output_dir"],sample=SAMPLES),
        perIndvcfGZ=expand("{path}/snp_calling/{sample}.vcf.gz",path=config["output_dir"],sample=SAMPLES)
    output:
        finalVCF=expand("{path}/snp_calling/snp.vcf.gz",path=config["output_dir"])
    params:
        outDir=expand("{path}/snp_calling/",path=config["output_dir"])
    conda:
        "src/env/bcftools.yaml"
    shell:
        """bcftools merge {params.outDir}/*.gz | bcftools view -i 'GT[*]="alt"&&REF!="N"' |  bgzip -c > {output.finalVCF}"""
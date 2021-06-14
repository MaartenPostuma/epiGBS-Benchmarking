configfile: "configRef.yaml"
import pandas as pd
import os
df = pd.read_csv(os.path.join(config["barcodes"]), sep='\t', dtype="object").set_index('Sample')
SAMPLES = df.index
FLAGS=["83","99","147","163"]


rule all:
    input:
        RTGsummary=expand("{out_dir}/RTG{score}/{sample}/summary.txt",sample=SAMPLES,out_dir=config["output_dir"],score=config["score_field"]),
        indAlignments=expand("{tmp}/{sample}-{flag}.bam",flag=FLAGS,tmp=config["tmp_dir"],sample=SAMPLES),
        test=expand("{tmp}/{sample}.depth",tmp=config["tmp_dir"],sample=SAMPLES),
        weighted_roc=expand("{out}/SquashRTG{score}/{sample}/weighted_roc.tsv.gz",out=config["output_dir"],score=config["score_field"],sample=SAMPLES),



        

rule sort:
    input:
        expand("{in_dir}/alignment/{{sample}}_trimmed_filt_merged.1_bismark_bt2_pe.bam",in_dir=config["input_dir"])
    output:
        bam=expand("{tmp}/{{sample}}_sorted.bam",tmp=config["tmp_dir"]),
    shell:
        """samtools sort {input} > {output.bam}"""
rule sort_index:
    input:
        bam=expand("{tmp}/{{sample}}_sorted.bam",tmp=config["tmp_dir"]),
    output:
        bai=expand("{tmp}/{{sample}}_sorted.bam.bai",tmp=config["tmp_dir"])
    shell:
        "samtools index {input.bam}"



rule extract_individual_alignments:
    input:
        bam=expand("{tmp}/{{sample}}_sorted.bam",tmp=config["tmp_dir"]),
        bai=expand("{tmp}/{{sample}}_sorted.bam.bai",tmp=config["tmp_dir"])
    output:
        temp(expand("{tmp}/{{sample}}-{{flag}}.bam",tmp=config["tmp_dir"]))
    params:
        flag="{flag}"
    shell:
        "samtools view -bf {params.flag} {input.bam} > {output}"

rule merge_alignments_CT:
    input:
        one=expand("{tmp}/{{sample}}-99.bam",tmp=config["tmp_dir"]),
        two=expand("{tmp}/{{sample}}-147.bam",tmp=config["tmp_dir"])
    output:
        temp(expand("{tmp}/{{sample}}.fixmate.CT.bam",tmp=config["tmp_dir"]))
    shell:
        """samtools merge - {input.one} {input.two} | samtools depth - -Q 10 | awk -v OFS="\t" '$3>0{{print $1,$2-1,$2,$3}}' > {output}"""

rule merge_alignments_GA:
    input:
        one=expand("{tmp}/{{sample}}-83.bam",tmp=config["tmp_dir"]),
        two=expand("{tmp}/{{sample}}-163.bam",tmp=config["tmp_dir"])
    output:
        temp(expand("{tmp}/{{sample}}.fixmate.GA.depth",tmp=config["tmp_dir"]))
    shell:
        """samtools merge - {input.one} {input.two} | samtools depth - -Q 10 | awk -v OFS="\t" '$3>0{{print $1,$2-1,$2,$3}}' > {output}"""

rule depth:
    input:
        GA=expand("{tmp}/{{sample}}.fixmate.GA.depth",tmp=config["tmp_dir"]),
        CT=expand("{tmp}/{{sample}}.fixmate.CT.bam",tmp=config["tmp_dir"])
    output:
        expand("{tmp}/{{sample}}.depth",tmp=config["tmp_dir"])
    shell:
        "bedtools intersect -a {input.CT} -b {input.GA} -u > {output}"


rule normRef:
    input:
        refIn="data/snp-calls/refAll.vcf.gz",
        ref=expand("{ref}",ref=config["ref"])
    output:
        ind=temp(expand("{tmp}/RefTemp.vcf.gz",tmp=config["tmp_dir"])),
        indTbi=temp(expand("{tmp}/RefTemp.vcf.gz.tbi",tmp=config["tmp_dir"]))
    shell:
        """
        bcftools view {input.refIn} | 
        bcftools norm -f {input.ref} -m- | 
        bcftools view -V indels,mnps,ref,bnd,other -Oz > {output.ind}
        bcftools tabix {output.ind} 
        """

rule indReference:
    input:
        ind=temp(expand("{tmp}/RefTemp.vcf.gz",tmp=config["tmp_dir"])),
        indTbi=temp(expand("{tmp}/RefTemp.vcf.gz.tbi",tmp=config["tmp_dir"]))
    output:
        ind=temp(expand("{tmp}/{{sample}}RefTemp.vcf.gz",tmp=config["tmp_dir"])),
        indTbi=temp(expand("{tmp}/{{sample}}RefTemp.vcf.gz.tbi",tmp=config["tmp_dir"]))
    params:
        sample='{sample}'
    shell:
        """
        bcftools view {input.ind} -s {params.sample} | bgzip -c > {output.ind}
        bcftools tabix {output.ind} 
        """



rule subset:
    input:
        depthIn=expand("{tmp}/{{sample}}.depth",tmp=config["tmp_dir"]),
        refIn=expand("{tmp}/{{sample}}RefTemp.vcf.gz",tmp=config["tmp_dir"])
    output:
        sampleSpecificRef=expand("{tmp}/{{sample}}Ref.vcf.gz",tmp=config["tmp_dir"])
    shell:
        """    
        cat <(bcftools view -h {input.refIn}) <(bedtools intersect -a {input.refIn} -b {input.depthIn} -u) | \
        bgzip -c > {output.sampleSpecificRef}
        """

rule tabixRef:
    input:
        sampleSpecificRef=expand("{tmp}/{{sample}}Ref.vcf.gz",tmp=config["tmp_dir"])
    output:
        sampleSpecificTabix=temp(expand("{tmp}/{{sample}}Ref.vcf.gz.tbi",tmp=config["tmp_dir"]))
    shell:
        "bcftools tabix {input.sampleSpecificRef} -f"

rule normalize:
    input:
        epiInput=expand("{in_dir}/snp_calling/snp.vcf.gz",in_dir=config["input_dir"]),
        ref=expand("{ref}",ref=config["ref"])
    output:
        norm=temp(expand("{tmp}/snp.norm.vcf.gz",tmp=config["tmp_dir"]))
    shell:
        "bcftools norm -f {input.ref} {input.epiInput} > {output.norm}"

rule filter:
    input:
        norm=expand("{tmp}/snp.norm.vcf.gz",tmp=config["tmp_dir"])
    output:
        filt=expand("{out}vcfs/{{sample}}.vcf.gz",out=config["output_dir"])
    params:
        sample="{sample}"
    shell:
        """cat <(bcftools view -h -s {params.sample} -V indels,mnps,ref,bnd,other {input.norm}) <(bcftools view -s {params.sample} -V indels,mnps,ref,bnd,other {input.norm} ) |  bgzip -c > {output.filt}"""

rule subset_calls:
    input:
        filt=expand("{out}vcfs/{{sample}}.vcf.gz",out=config["output_dir"]),
        depthIn=expand("{tmp}/{{sample}}.depth",tmp=config["tmp_dir"])
    output:
        filt=expand("{tmp}/{{sample}}.vcf.gz",tmp=config["output_dir"])
    shell:
        """
        cat <(bcftools view -h {input.filt}) <(bedtools intersect -a {input.filt} -b {input.depthIn} -u) | \
        bgzip -c > {output.filt}
        
        """
rule tabix_filter:
    input:
        filt=expand("{tmp}/{{sample}}.vcf.gz",tmp=config["output_dir"])
    output:
        filt=expand("{tmp}/{{sample}}.vcf.gz.tbi",tmp=config["output_dir"])
    shell:
        "bcftools tabix {input.filt}"

rule referenceFormat:
    input:
        ref=expand("{ref}",ref=config["ref"])
    output:
        refSDF=expand("{ref}.sdf",ref=config["ref"])
    shell:
        "rtg format -o {output.refSDF} {input.ref}"

rule high_confidence:
    input:
        sampleSpecificRef=expand("{tmp}/{{sample}}Ref.vcf.gz",tmp=config["tmp_dir"]),
        sampleSpecificTabix=expand("{tmp}/{{sample}}Ref.vcf.gz.tbi",tmp=config["tmp_dir"])
    output:
        high_confidence=expand("{tmp}/{{sample}}_high_confidence.bed",tmp=config["tmp_dir"])
    shell:
        """bcftools view -e 'GT=="./." || GT=="./1" || GT=="0/."' {input.sampleSpecificRef} | bedtools merge -i stdin > {output.high_confidence}"""

rule rtg:
    input:
        sampleSpecificRef=expand("{tmp}/{{sample}}Ref.vcf.gz",tmp=config["tmp_dir"]),
        sampleSpecificTabix=expand("{tmp}/{{sample}}Ref.vcf.gz.tbi",tmp=config["tmp_dir"]),
        filt=expand("{tmp}/{{sample}}.vcf.gz",tmp=config["output_dir"]),
        tabix=expand("{tmp}/{{sample}}.vcf.gz.tbi",tmp=config["output_dir"]),
        refSDF=expand("{ref}.sdf",ref=config["ref"]),
        high_confidence=expand("{tmp}/{{sample}}_high_confidence.bed",tmp=config["tmp_dir"])
    params:
        outDirTemp=expand("{out}/RTG/{{sample}}temp/",out=config["output_dir"]),
        outDir=expand("{out}/RTG{score}/{{sample}}/",out=config["output_dir"],score=config["score_field"]),
        sample="{sample}",
        vcfScoreField=config["score_field"]
    output:
        summary=expand("{out}/RTG{score}/{{sample}}/summary.txt",out=config["output_dir"],score=config["score_field"]),
        fp=expand("{out}/RTG{score}/{{sample}}/fp.vcf.gz",out=config["output_dir"],score=config["score_field"]),
        fn=expand("{out}/RTG{score}/{{sample}}/fn.vcf.gz",out=config["output_dir"],score=config["score_field"]),
        tp=expand("{out}/RTG{score}/{{sample}}/tp.vcf.gz",out=config["output_dir"],score=config["score_field"]),
        weighted_roc=expand("{out}/RTG{score}/{{sample}}/weighted_roc.tsv.gz",out=config["output_dir"],score=config["score_field"]),
    shell:
        """
        mkdir {params.outDir} -p
        rtg vcfeval -b {input.sampleSpecificRef} -c {input.filt} -t {input.refSDF} -o {params.outDirTemp} \
        --sample={params.sample} \
        --vcf-score-field={params.vcfScoreField} \
        --evaluation-regions={input.high_confidence} 
        mv {params.outDirTemp}/* {params.outDir}
        rm -r {params.outDirTemp}
        """


rule rtg_squash:
    input:
        sampleSpecificRef=expand("{tmp}/{{sample}}Ref.vcf.gz",tmp=config["tmp_dir"]),
        sampleSpecificTabix=expand("{tmp}/{{sample}}Ref.vcf.gz.tbi",tmp=config["tmp_dir"]),
        filt=expand("{tmp}/{{sample}}.vcf.gz",tmp=config["output_dir"]),
        tabix=expand("{tmp}/{{sample}}.vcf.gz.tbi",tmp=config["output_dir"]),
        refSDF=expand("{ref}.sdf",ref=config["ref"]),
        high_confidence=expand("{tmp}/{{sample}}_high_confidence.bed",tmp=config["tmp_dir"])
    params:
        outDirTemp=expand("{out}/SquashRTG/{{sample}}temp/",out=config["output_dir"]),
        outDir=expand("{out}/SquashRTG{score}/{{sample}}/",out=config["output_dir"],score=config["score_field"]),
        sample="{sample}",
        vcfScoreField=config["score_field"]
    output:
        summary=expand("{out}/SquashRTG{score}/{{sample}}/summary.txt",out=config["output_dir"],score=config["score_field"]),
        fp=expand("{out}/SquashRTG{score}/{{sample}}/fp.vcf.gz",out=config["output_dir"],score=config["score_field"]),
        fn=expand("{out}/SquashRTG{score}/{{sample}}/fn.vcf.gz",out=config["output_dir"],score=config["score_field"]),
        tp=expand("{out}/SquashRTG{score}/{{sample}}/tp.vcf.gz",out=config["output_dir"],score=config["score_field"]),
        weighted_roc=expand("{out}/SquashRTG{score}/{{sample}}/weighted_roc.tsv.gz",out=config["output_dir"],score=config["score_field"]),
    shell:
        """
        mkdir {params.outDir} -p
        rtg vcfeval -b {input.sampleSpecificRef} -c {input.filt} -t {input.refSDF} -o {params.outDirTemp} \
        --sample={params.sample} \
        --vcf-score-field={params.vcfScoreField} \
        --evaluation-regions={input.high_confidence} \
        --squash-ploidy
        mv {params.outDirTemp}/* {params.outDir}
        rm -r {params.outDirTemp}
        """
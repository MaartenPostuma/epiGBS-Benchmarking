configfile: "configRef.yaml"
import pandas as pd
import os
df = pd.read_csv(os.path.join(config["barcodes"]), sep='\t', dtype="object").set_index('Sample')
SAMPLES = df.index

rule all:
    input:
        RTGsummary=expand("{out}/RTG{score}/{sample}/summary.txt",sample=SAMPLES,out=config["output_dir"],score=config["score_field"]),
        

rule sort:
    input:
        expand("{in_dir}/alignment/{{sample}}_trimmed_filt_merged.1_bismark_bt2_pe.bam",in_dir=config["input_dir"])
    output:
        temp(expand("{tmp}/{{sample}}_sorted.bam",tmp=config["tmp_dir"]))
    shell:
        "samtools sort {input} > {output}"

rule depth:
    input:
        expand("{tmp}/{{sample}}_sorted.bam",tmp=config["tmp_dir"])
    output:
        temp(expand("{tmp}/{{sample}}.depth",tmp=config["tmp_dir"]))
    shell:
        "samtools depth {input} |  awk '$3>0{{print $1,$2-1,$2,$3}}' | sed 's/ /\t/g' > {output}"

rule subset:
    input:
        depthIn=expand("{tmp}/{{sample}}.depth",tmp=config["tmp_dir"]),
        refIn="data/snp-calls/refAll.vcf.gz"
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
        filt=expand("{out}/vcfs/{{sample}}.vcf.gz",out=config["output_dir"])
    params:
        sample="{sample}"
    shell:
        """cat <(bcftools view -h -s {params.sample} -V indels,mnps,ref,bnd,other {input.norm}) <(bcftools view -s {params.sample} -V indels,mnps,ref,bnd,other {input.norm}  | awk '$0~"^#" || ($10~"^1/1" || $10~"^1/0" || $10~"^0/1")' )| bgzip -c > {output.filt}"""

rule tabix_filter:
    input:
        filt=expand("{out}/vcfs/{{sample}}.vcf.gz",out=config["output_dir"])
    output:
        filt=expand("{out}/vcfs/{{sample}}.vcf.gz.tbi",out=config["output_dir"])
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
        filt=expand("{out}/vcfs/{{sample}}.vcf.gz",out=config["output_dir"]),
        tabix=expand("{out}/vcfs/{{sample}}.vcf.gz.tbi",out=config["output_dir"]),
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
        filt=expand("{out}/vcfs/{{sample}}.vcf.gz",out=config["output_dir"]),
        tabix=expand("{out}/vcfs/{{sample}}.vcf.gz.tbi",out=config["output_dir"]),
        refSDF=expand("{ref}.sdf",ref=config["ref"])
    params:
        outDirTemp=expand("{out}/RTGSquash/{{sample}}temp/",out=config["output_dir"]),
        outDir=expand("{out}/RTGSquash{score}/{{sample}}/",out=config["output_dir"],score=config["score_field"]),
        sample="{sample}",
        vcfScoreField=config["score_field"]
    output:
        summary=expand("{out}/RTGSquash{score}/{{sample}}/summary.txt",out=config["output_dir"],score=config["score_field"]),
        fp=expand("{out}/RTGSquash{score}/{{sample}}/fp.vcf.gz",out=config["output_dir"],score=config["score_field"]),
        fn=expand("{out}/RTGSquash{score}/{{sample}}/fn.vcf.gz",out=config["output_dir"],score=config["score_field"]),
        tp=expand("{out}/RTGSquash{score}/{{sample}}/tp.vcf.gz",out=config["output_dir"],score=config["score_field"]),
        weighted_roc=expand("{out}/RTGSquash{score}/{{sample}}/weighted_roc.tsv.gz",out=config["output_dir"],score=config["score_field"]),
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
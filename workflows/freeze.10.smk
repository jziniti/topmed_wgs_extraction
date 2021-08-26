CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

FREEZE10 = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10a/minDP0/freeze.10a.chr{chrom}.pass_and_fail.gtonly.minDP0.bcf'
FREEZE10_IRC = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10b.irc/minDP0/freeze.10a.chr{chrom}.pass_and_fail.gtonly.minDP0'
LOWQUAL = '/proj/regeps/regep00/studies/CAMP/data/dna/whole_genome/TopMed/data/freezes/lowqual/'

STUDIES = config['studies']
TARGETS = expand('tmp/{s_studyid}.{chr}.extract.bcf', s_studyid=STUDIES, chr=CHROMOSOMES)

rule done: input: TARGETS

rule stashq:
    input: "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10a/manifests/nwdids.txt"
    output: temp("tmp/all.stashq.txt")
    conda: "../envs/stashq.yaml"
    shell: "stashq alias `cat {input}` > {output}"

rule get_study_ids:
    input: "tmp/all.stashq.txt"
    output: nwdids=temp("tmp/{s_studyid}.nwds.txt")
    conda: "../envs/sapphire8.yaml"
    shell: "grep {wildcards.s_studyid} {input} | cut -d \" \" -f 1 > {output}"

rule sapphire_nwd_list:
    input: 
    output: nwdids=temp("tmp/{s_studyid}.sapphire.txt")
    conda: "../envs/sapphire8.yaml"
    params:
        alias_prefix='NWD',
        s_studyid='{s_studyid}',
    script: "../scripts/python/read_sapphire_aliases.py"

rule extract:
    input:
        samples="tmp/{s_studyid}.nwds.txt",
        bcf="/proj/edith/regeps/regep00/studies/COPDGene/data/wgs/TopMed/data/freezes/freeze.10a/minDP10/freeze.10a.chr{chr}.pass_and_fail.gtonly.minDP10.bcf",
    output:
        bcf=temp("tmp/{s_studyid}.{chr}.PASS.bcf")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools view -S {input.samples} -i 'FILTER=\"PASS\"' -c 1 -O b -m2 -M2 --force-samples --types snps {input.bcf} > {output.bcf}"

rule setid:
    input:
        stashq="tmp/{s_studyid}.stashq",
        vcf=rules.extract.output.bcf,
    output: temp("tmp/{s_studyid}.{chr}.annotated.bcf")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O b -o {output} {input}"

rule convert_to_plink:
    input: bcf=rules.setid.output.bcf
    output:
        bed=temp("tmp/{s_studyid}_annotated_plink_chr{chr}.bed"),
        bim=temp("tmp/{s_studyid}_annotated_plink_chr{chr}.bim"),
        fam=temp("tmp/{s_studyid}_annotated_plink_chr{chr}.fam"),
    conda: "../envs/plink.yaml"
    shell:
        """plink --bcf {input} --make-bed --out tmp/{wildcards.s_studyid}_annotated_plink_chr{wildcards.chrom}"""
        
rule merge:
    input:
        expand("{s_studyid}_annotated_plink_chr{chr}.fam", chr=CHROMOSOMES),
        expand("{s_studyid}_annotated_plink_chr{chr}.bim", chr=CHROMOSOMES),
        expand("{s_studyid}_annotated_plink_chr{chr}.bed", chr=CHROMOSOMES),
        merge_list='tmp/{s_studyid}.merge_list',
    output:
        bed="tmp/{s_studyid}_annotated_plink_merged.bed",
        bim="tmp/{s_studyid}_annotated_plink_merged.bim",
        fam="tmp/{s_studyid}_annotated_plink_merged.fam"
    conda: "../envs/plink.yaml"
    params:
        geno=float(config.get('geno_filter', 0.02)),
        maf=float(config.get('maf_filter', 0.0001)),
    shell: """plink --merge-list {input.merge_list} --geno {params.geno} --maf {params.maf} --make-bed --out tmp/{wildcards.s_studyid}_annotated_plink_merged"""

rule get_pedigree:
    input: fam=rules.merge.output.fam
    output: "tmp/{s_studyid}_annotated_plink_merged_ped.fam"
    conda: "../envs/r.yaml"
    script: "../scripts/R/put_in_pedigree.R"
    
rule sex_check:
    input:
        bed=rules.merge.output.bed,
        bim=rules.merge.output.bim,
        fam=rules.merge.output.fam,
    output: temp("tmp/{s_studyid}.sexcheck")
    conda: "../envs/plink.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --check-sex --out tmp/{wildcards.s_studyid}"

rule king:
    input:
        bed=rules.merge.output.bed,
        bim=rules.merge.output.bim,
        fam=rules.merge.output.fam,
    output: temp("tmp/{s_studyid}_king_duplicate.con")
    conda: "../envs/king.yaml"
    shell: "king -b {input.bed} --fam {input.fam} --bim {input.bim} --duplicate --prefix tmp/{wildcards.s_studyid}_king_duplicate"

#rule variant_qc:
#    input: "tmp/{s_studyid}.annotate.vcf.gz"
##    output: temp("tmp/{s_studyid}.variants.keep")
#    conda: "../envs/bcftools.yaml"
#    shell: "bcftools"
#
#rule sample_qc:
#    input: "tmp/{s_studyid}.annotate.vcf.gz"
#    output: temp("tmp/{s_studyid}.samples.keep")
#    conda: "../envs/bcftools.yaml"
#    shell: "bcftools"
#
#rule filter:
#    input:
#        "tmp/{s_studyid}.annotate.vcf.gz",
#        "tmp/{s_studyid}.samples.keep",
#        "tmp/{s_studyid}.variants.keep",
#    output: temp("{s_studyid}.bed")
#    conda: "../envs/bcftools.yaml"
#    shell: "bcftools"

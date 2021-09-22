CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

FREEZE9B = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.9b/minDP10/freeze.9b.chr{chrom}.pass_and_fail.gtonly.minDP10.bcf'
FREEZE10 = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10a/minDP10/freeze.10a.chr{chrom}.pass_and_fail.gtonly.minDP10.bcf'
FREEZE10_IRC = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10b.irc/minDP10/freeze.10b.chr{chrom}.pass_and_fail.gtonly.minDP10.bcf'

TMP = Path('tmp')

STUDIES = sorted(config.keys())
TARGETS = expand('tmp/{s_studyid}.{chrom}.PASS.bcf', s_studyid=STUDIES, chrom=CHROMOSOMES)
# TARGETS += expand('tmp/{s_studyid}_king.kin0', s_studyid=STUDIES)
TARGETS += expand("tmp/{s_studyid}.sexcheck", s_studyid=STUDIES)
TARGETS += expand("tmp/{s_studyid}_reference.stashq.txt", s_studyid=STUDIES)
TARGETS += expand("tmp/{s_studyid}_king_duplicate.con", s_studyid=STUDIES)
TARGETS += expand(TMP/'{s_studyid}_annotated_plink_merged.{suffix}', s_studyid=STUDIES, suffix=('frq', 'imiss', 'lmiss', 'het', 'hwe', 'sexcheck'))

BCF_PATH_FOR_STUDY = {
    'GECOPD': FREEZE10,
    'EOCOPD': FREEZE10_IRC,
    'ECLPSE': FREEZE10_IRC,
    'GLAXO':  FREEZE10_IRC,
    'PLCOPD': FREEZE10_IRC,
    'LTCOPD': FREEZE10_IRC,
    'LTRC': FREEZE10_IRC,
}

SOURCE_FREEZE_FOR_STUDY = {
    ## 'GECOPD': '10a.irc',
    'EOCOPD': '10a.irc',
    'ECLPSE': '10a.irc',
    'GLAXO':  '10a.irc',
    'PLCOPD': '10a.irc',
    'LTCOPD': '10a.irc',
    'LTRC': '10a.irc',
}

include: "gwas_qc_pipeline.smk"
# include: "methylation_ref_concordance.smk"
# include: "test_notebooks.smk"

rule done: input: TARGETS

rule stashq:
    input: "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.{freeze}/manifests/nwdids.txt"
    output: temp("tmp/wgs.{freeze}.stashq.txt")
    conda: "../envs/stashq.yaml"
    shell: "stashq alias `cat {input}` > {output}"

rule stashq_reference:
    input: lambda w: f"{config[w.s_studyid]['known_good_reference']}.fam"
    output: temp("tmp/{s_studyid}_reference.stashq.txt")
    conda: "../envs/stashq.yaml"
    shell: "stashq alias `awk '{{print $2}}' {input}` > {output}"

rule get_study_ids:
    input: lambda w: f"tmp/wgs.{SOURCE_FREEZE_FOR_STUDY.get(w.s_studyid, '10a')}.stashq.txt"
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
        #bcf=lambda w: f"/proj/edith/regeps/regep00/studies/COPDGene/data/wgs/TopMed/data/freezes/freeze.{SOURCE_FREEZE_FOR_STUDY.get(w.s_studyid, '10a')}/minDP10/freeze.10a.chr{{chrom}}.pass_and_fail.gtonly.minDP10.bcf",
        bcf=lambda w: BCF_PATH_FOR_STUDY.get(w.s_studyid, FREEZE10),
    output:
        bcf=temp("tmp/{s_studyid}.{chrom}.PASS.bcf")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools view -S {input.samples} -i 'FILTER=\"PASS\"' -c 1 -O b -m2 -M2 --force-samples --types snps {input.bcf} -o {output.bcf}"

#rule plink_freq:
#    input: bcf="tmp/{s_studyid}.{chrom}.PASS.bcf"
#    output: frq=temp("tmp/{s_studyid}.{chrom}.PASS.afreq")
#    conda: "../envs/plink2a.yaml"
#    shell: "plink2 --bcf {input.bcf} --out tmp/{wildcards.s_studyid}.{wildcards.chrom}.PASS --freq"
    
rule setid:
    input:
        vcf=rules.extract.output.bcf,
    output: bcf=temp("tmp/{s_studyid}.{chrom}.annotated.bcf")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O b -o {output} {input.vcf}"

rule convert_to_plink:
    input: bcf=rules.setid.output.bcf
    output:
        bed=temp("tmp/{s_studyid}_annotated_plink_chr{chrom}.bed"),
        bim=temp("tmp/{s_studyid}_annotated_plink_chr{chrom}.bim"),
        fam=temp("tmp/{s_studyid}_annotated_plink_chr{chrom}.fam"),
    conda: "../envs/plink2a.yaml"
    shell:
        """plink --allow-extra-chr  --bcf {input} --make-bed --out tmp/{wildcards.s_studyid}_annotated_plink_chr{wildcards.chrom}"""
        
rule merge_list:
    input:
    output: merge_list=temp('tmp/{s_studyid}.merge_list')
    run:
        with open(f'tmp/{wildcards.s_studyid}.merge_list', 'w') as fh:
            for chrom in CHROMOSOMES:
                fh.write(f'tmp/{wildcards.s_studyid}_annotated_plink_chr{chrom}.bed tmp/{wildcards.s_studyid}_annotated_plink_chr{chrom}.bim tmp/{wildcards.s_studyid}_annotated_plink_chr{chrom}.fam\n')

rule merge:
    input:
        expand("tmp/{{s_studyid}}_annotated_plink_chr{chrom}.fam", chrom=CHROMOSOMES),
        expand("tmp/{{s_studyid}}_annotated_plink_chr{chrom}.bim", chrom=CHROMOSOMES),
        expand("tmp/{{s_studyid}}_annotated_plink_chr{chrom}.bed", chrom=CHROMOSOMES),
        merge_list='tmp/{s_studyid}.merge_list',
    output:
        bed="tmp/{s_studyid}_annotated_plink_merged.bed",
        bim="tmp/{s_studyid}_annotated_plink_merged.bim",
        fam="tmp/{s_studyid}_annotated_plink_merged.fam"
    conda: "../envs/plink2a.yaml"
    params:
        geno=float(config.get('geno_filter', 0.02)),
        maf=float(config.get('maf_filter', 0.0001)),
    shell: """plink --pmerge-list {input.merge_list} --geno {params.geno} --maf {params.maf} --make-bed --out tmp/{wildcards.s_studyid}_annotated_plink_merged"""

rule get_pedigree:
    input: fam=rules.merge.output.fam
    output: "tmp/{s_studyid}_annotated_plink_merged_ped.fam"
    conda: "../envs/r.yaml"
    script: "../scripts/R/put_in_pedigree.R"
    
#rule sex_check:
#    input:
#        bed=rules.merge.output.bed,
#        bim=rules.merge.output.bim,
#        fam=rules.merge.output.fam,
#    output: temp("tmp/{s_studyid}.sexcheck")
#    conda: "../envs/plink.yaml"
#    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --check-sex --out tmp/{wildcards.s_studyid}"

rule king:
    input:
        bed=rules.merge.output.bed,
        bim=rules.merge.output.bim,
        fam=rules.merge.output.fam,
    output: temp("tmp/{s_studyid}_king_duplicate.con")
    conda: "../envs/king.yaml"
    shell: "king -b {input.bed} --fam {input.fam} --bim {input.bim} --duplicate --prefix tmp/{wildcards.s_studyid}_king_duplicate"

rule concordance_vs_reference:
    input:  
        qbed=rules.merge.output.bed,
        qbim=rules.merge.output.bim,
        qfam=rules.merge.output.fam,
        rbed=lambda w: f"{config[w.s_studyid]['known_good_reference']}.bed",
        rbim=lambda w: f"{config[w.s_studyid]['known_good_reference']}.bim",
        rfam=lambda w: f"{config[w.s_studyid]['known_good_reference']}.fam",
    output: temp("tmp/{s_studyid}_king.kin0"),
    conda: "../envs/king.yaml"
    params:
        Q="tmp/{s_studyid}_king",
        R=lambda w: "{config[w.s_studyid]['known_good_reference']}",
    shell: "king -b {params.Q},{params.R} --prefix {params.Q} --kinship"

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

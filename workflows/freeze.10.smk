CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

FREEZE9B = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.9b/minDP10/freeze.9b.chr{chrom}.pass_and_fail.gtonly.minDP10.bcf'
FREEZE10 = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10a/minDP10/freeze.10a.chr{chrom}.pass_and_fail.gtonly.minDP10.bcf'
FREEZE10_IRC = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10b.irc/minDP10/freeze.10b.chr{chrom}.pass_and_fail.gtonly.minDP10.bcf'

TMP = Path('tmp')
NOTEBOOKS = Path('notebook_logs')
FLAGS = Path('flags')

STUDIES = sorted(config.keys())
TARGETS = []
for s_studyid in STUDIES:
    if config[s_studyid].get('run_reference_concordance', False):
        TARGETS.append(f'tmp/{s_studyid}_reference_concordance.con')
    if config[s_studyid].get('run_methylation_concordance', False):
        TARGETS.append(TMP/f'{s_studyid}_methylation_freeze10.kin0')
    if config[s_studyid].get('run_qc_notebook', False):
        TARGETS.append(FLAGS/f'{s_studyid}_samples.csv')
        TARGETS.append(FLAGS/f'{s_studyid}_markers.csv')
    #if config[s_studyid].get('run_pca_pipeline', False):
    #    TARGETS.append(TMP/f'{s_studyid}_ready_for_umich.done')
    #if config[s_studyid].get('run_rna_concordance', False):
    #    TARGETS.append(TMP/f'{s_studyid}_rna_king_results_summary.html')
    #    TARGETS.append(TMP/f'{s_studyid}_rna_king_results_summary.csv')
    TARGETS.append(f'multiomics/{s_studyid}/ANNOTATED_MANIFEST.csv')
# TARGETS += expand("tmp/{s_studyid}.sexcheck", s_studyid=STUDIES)
# TARGETS += expand("tmp/{s_studyid}_reference.stashq.txt", s_studyid=STUDIES)
# TARGETS += expand("tmp/{s_studyid}_king_duplicate.con", s_studyid=STUDIES)
# TARGETS += expand(TMP/'{s_studyid}_annotated_plink_merged.{suffix}', s_studyid=STUDIES, suffix=('frq', 'imiss', 'lmiss', 'het', 'hwe', 'sexcheck'))


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
    'GECOPD': '10a',
    'EOCOPD': '10a.irc',
    'ECLPSE': '10a.irc',
    'GLAXO':  '10a.irc',
    'PLCOPD': '10a.irc',
    'LTCOPD': '10a.irc',
    'LTRC': '10a.irc',
}

### Modules we will be using
include: "gwas_qc_pipeline.smk"
include: "methylation_ref_concordance.smk"
include: "reaka_test_notebooks.smk"
include: "multiomics.smk"
include: "plate103.smk"
#include: "rnaseq_concordance.smk"
include: "cdnm_pca_pipeline_shim.smk"
include: "duplicate_vs_reference.smk"

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

rule setid:
    input:
        vcf=rules.extract.output.bcf,
    output: bcf=temp("tmp/{s_studyid}.{chrom}.annotated.bcf")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -O b -o {output} {input.vcf}"

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

rule merge_and_filter:
    input:
        expand("tmp/{{s_studyid}}_annotated_plink_chr{chrom}.fam", chrom=CHROMOSOMES),
        expand("tmp/{{s_studyid}}_annotated_plink_chr{chrom}.bim", chrom=CHROMOSOMES),
        expand("tmp/{{s_studyid}}_annotated_plink_chr{chrom}.bed", chrom=CHROMOSOMES),
        merge_list='tmp/{s_studyid}.merge_list',
    output:
        bed="tmp/{s_studyid}_annotated_plink_merged_noPAR.bed",
        bim="tmp/{s_studyid}_annotated_plink_merged_noPAR.bim",
        fam="tmp/{s_studyid}_annotated_plink_merged_noPAR.fam"
    conda: "../envs/plink2a.yaml"
    params:
        geno=float(config.get('geno_filter', 0.02)),
        maf=float(config.get('maf_filter', 0.0001)),
    shell: """plink --pmerge-list {input.merge_list} --geno {params.geno} --maf {params.maf} --make-bed --out tmp/{wildcards.s_studyid}_annotated_plink_merged_noPAR"""

rule make_PAR_list:
    input:
        bim="tmp/{s_studyid}_annotated_plink_merged_noPAR.bim",
    output:
        TMP/"{s_studyid}_par_newchr.txt"
    params:
        par1_start=10001,
        par1_end=2781479,
        par2_start=155701383,
        par2_end=156030895,
        chrX_label="X"
    shell:
        "awk '($1 == \"{params.chrX_label}\") && (($4 >= {params.par1_start} && $4 <= {params.par1_end}) || ($4 >= {params.par2_start} && $4 <= {params.par2_end})) {{print $2\"\\t25\"}}' {input} > {output}"

rule annotate_PAR:
    input:
        bed=TMP/"{s_studyid}_annotated_plink_merged_noPAR.bed",
        bim=TMP/"{s_studyid}_annotated_plink_merged_noPAR.bim",
        fam=TMP/"{s_studyid}_annotated_plink_merged_noPAR.fam",
        parlist=TMP/"{s_studyid}_par_newchr.txt",
    output:
        bed=TMP/"{s_studyid}_annotated_plink_merged.bed",
        bim=TMP/"{s_studyid}_annotated_plink_merged.bim",
        fam=TMP/"{s_studyid}_annotated_plink_merged.fam"
    conda: "../envs/plink.yaml"
    shell: """plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --update-chr {input.parlist} --out tmp/{wildcards.s_studyid}_annotated_plink_merged --make-bed"""

rule convert_to_vcf:
    input:
        bed="tmp/{s_studyid}_annotated_plink_merged.bed",
        bim="tmp/{s_studyid}_annotated_plink_merged.bim",
        fam="tmp/{s_studyid}_annotated_plink_merged.fam"
    output:
        bed="tmp/{s_studyid}_annotated_plink_merged.vcf.gz",
    conda: "../envs/plink2a.yaml"
    shell: """plink --bfile tmp/{wildcards.s_studyid}_annotated_plink_merged --out tmp/{wildcards.s_studyid}_annotated_plink_merged --recode vcf bgz"""

rule get_pedigree:
    input: fam=rules.merge_and_filter.output.fam
    output: "tmp/{s_studyid}_annotated_plink_merged_ped.fam"
    conda: "../envs/r.yaml"
    script: "../scripts/R/put_in_pedigree.R"
    
rule king_duplicate_check:
    input:
        bed=rules.merge_and_filter.output.bed,
        bim=rules.merge_and_filter.output.bim,
        fam=rules.merge_and_filter.output.fam,
    output: temp("tmp/{s_studyid}_king_duplicate.con")
    conda: "../envs/king.yaml"
    shell: "king -b {input.bed} --duplicate --prefix tmp/{wildcards.s_studyid}_king_duplicate"

rule filter:
    input:
        bed="tmp/{s_studyid}.annotated_plink_merged.bed",
        bim="tmp/{s_studyid}.annotated_plink_merged.bim",
        fam="tmp/{s_studyid}.annotated_plink_merged.fam",
        rem="flags/{s_studyid}_samples.csv",
        exclude="flages/{s_studyid}_markers.csv",
    output:
        bed="output_freezes/{s_study}/{s_studyid}_freeze.10.bed",
        bim="output_freezes/{s_study}/{s_studyid}_freeze.10.bim",
        fam="output_freezes/{s_study}/{s_studyid}_freeze.10.fam",
    conda: "../envs/plink2a.yaml"
    shell: "plink --bfile  tmp/{wildcards.s_studyid}.annotated_plink_merged.bed \
                  --out output_freezes/{wildcards.s_study}/{wildcards.s_studyid}_freeze.10 \
                  --remove {input.rem} \
                  --exclude {input.exclude} \
                  --make-bed"

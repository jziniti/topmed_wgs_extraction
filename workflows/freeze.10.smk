import datetime
from snakemake.utils import validate
import sys

## validate(config, "../schema/config.schema.yaml")

CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

FREEZE9B = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.9b/minDP10/freeze.9b.chr{chrom}.pass_and_fail.gtonly.minDP10.bcf'
FREEZE10 = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10a/minDP10/freeze.10a.chr{chrom}.pass_and_fail.gtonly.minDP10.bcf'
FREEZE10_IRC = '/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10b.irc/minDP10/freeze.10b.chr{chrom}.pass_and_fail.gtonly.minDP10.bcf'

CMDLINE = ' '.join(sys.argv)
TIMESTAMP = datetime.datetime.now().strftime('%Y%m%d%H%M%S')

TMP = Path('tmp')
NOTEBOOKS = Path('notebook_logs')
FLAGS = Path('flags')

STUDIES = sorted(config.keys())
TARGETS = []
for s_studyid in STUDIES:
    # TARGETS.append(f"output_freezes/{s_studyid}/{s_studyid}_freeze.10.bed")
    if config[s_studyid].get('run_multiomics_qc', False):
        TARGETS.append(f'multiomics/{s_studyid}/ANNOTATED_MANIFEST.csv')
    if config[s_studyid].get('run_reference_concordance', False):
        TARGETS.append(f'tmp/{s_studyid}_reference_concordance.con')
    if config[s_studyid].get('run_methylation_concordance', False):
        TARGETS.append(TMP/f'{s_studyid}_methylation_freeze10.kin0')
    #if config[s_studyid].get('run_pca_pipeline', False):
    #    TARGETS.append(TMP/f'{s_studyid}_ready_for_umich.done')
    #if config[s_studyid].get('run_rna_concordance', False):
    #    TARGETS.append(TMP/f'{s_studyid}_rna_king_results_summary.html')
    #    TARGETS.append(TMP/f'{s_studyid}_rna_king_results_summary.csv')
    TARGETS.append(f'output_freezes/{s_studyid}/{s_studyid}_freeze.10.pep.yaml')
    TARGETS.append(f'output_freezes/{s_studyid}/{s_studyid}_freeze.10.pep.csv')
    TARGETS.append(f'output_freezes/{s_studyid}/md5sums.txt')
    for chrom in CHROMOSOMES:
        TARGETS.append(f'output_freezes/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.bcf')
        TARGETS.append(f'output_freezes/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.bcf.csi')

## I'm not sure if these are needed anymore ...
#TARGETS.append(TMP/'IRCALL_exome_duplicate.con')
#TARGETS.append(TMP/'IRCALL_axiom1_duplicate.con') ### This one segfaults ... datasets too big?
#TARGETS.append(TMP/'IRCALL_GECOPD_duplicate.con')
#TARGETS.append(TMP/'IRCALL_CAMP_duplicate.con')
#TARGETS.append(TMP/'IRCALL_CRA_duplicate.con')
#TARGETS.append(TMP/'IRCALL_ECLPSE_duplicate.con')
#TARGETS.append(TMP/"COPD_lost_samples_rendered.Rmd")

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

wildcard_constraints:
    s_studyid='[A-Z]+'

### Modules we will be using
if "CAMP" in config.keys():
    include: "camp.smk"

if 'GECOPD' in config.keys():
    include: "gecopd_heterozygosity.smk"
    include: "plate103.smk"
    # include: "ircall_concordance.smk"
    # include: "../COPD_lost_samples/Snakefile"
    # include: "gecopd_all_rna.smk"
    # include: "freeze3_dbgap_submission.smk"
    ## include: "rnaseq_concordance.smk"
if 'ECLPSE' in config.keys():
    include: "eclipse.smk"

include: "gwas_qc_pipeline.smk"
include: "methylation_ref_concordance.smk"
include: "run_qc_notebook.smk"
include: "multiomics.smk"
include: "duplicate_vs_reference.smk"
include: "eocopd_exome_sexcheck.smk"
include: "cra_forced_good_sexcheck.smk"
include: "glaxo_reference_prep.smk"
include: "exome6k_preparation.smk"

#include: "cdnm_pca_pipeline_shim.smk"

rule all: input: TARGETS

rule stashq:
    input: "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.{freeze}/manifests/nwdids.txt"
    output: temp("tmp/wgs.{freeze}.stashq.txt")
    conda: "../envs/stashq.yaml"
    ## shell: "stashq alias `cat {input}` > {output}"
    params: script=srcdir("../scripts/bash/split_stashq.sh")
    shell: "bash {params.script} {input} > {output}"

rule stashq_reference:
    input: lambda w: f"{config[w.s_studyid]['known_good_reference']}.fam"
    output: temp("tmp/{s_studyid}_reference.stashq.txt")
    conda: "../envs/stashq.yaml"
    shell: "stashq alias `awk '{{print $2}}' {input}` > {output}"

# ruleorder: get_aliases_for_GECOPD > get_aliases_for_study
# rule get_aliases_for_GECOPD:
#    input:
#        stashq=lambda w: f"tmp/wgs.{SOURCE_FREEZE_FOR_STUDY.get(w.s_studyid, '10a')}.stashq.txt",
#        recovered_samples="code/COPD_lost_samples/output/assigned.csv",
#        recovered_sk38="code/COPD_lost_samples/output/assigned_from_sk38.txt",
#    output: nwdids=temp("tmp/{s_studyid}.nwds.txt")
#    conda: "../envs/sapphire8.yaml"
#    shell: """grep {wildcards.s_studyid} {input.stashq} | cut -d \" \" -f 1 > {output}; \
#            awk -F , '$2 == \"{wildcards.s_studyid}\" {{ print $1}}' {input.recovered_samples} >> {output}; \
#            grep -v NWDID {input.recovered_sk38} >> {output}"""

rule get_aliases_for_study:
    input:
        stashq=lambda w: f"tmp/wgs.{SOURCE_FREEZE_FOR_STUDY.get(w.s_studyid, '10a')}.stashq.txt",
        recovered_samples="code/COPD_lost_samples/output/assigned.csv",
    output: nwdids=temp("tmp/{s_studyid}.nwds.txt")
    conda: "../envs/sapphire8.yaml"
    shell: "grep {wildcards.s_studyid} {input.stashq} | cut -d \" \" -f 1 > {output}; awk -F , '$2 == \"{wildcards.s_studyid}\" {{ print $1}}' {input.recovered_samples} >> {output}"

rule freeze10_irc_manifest:
    input: "tmp/wgs.10a.irc.stashq.txt"
    output: nwdids=temp("tmp/IRCALL.nwds.txt")
    conda: "../envs/sapphire8.yaml"
    shell: "cut -d \" \" -f 1 {input} > {output}"

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
        #bcf=protected("tmp/{s_studyid}.{chrom}.PASS.bcf")
        bcf="tmp/{s_studyid}.{chrom}.PASS.bcf"
    conda: "../envs/bcftools.yaml"
    shell: "bcftools view -S {input.samples} -i 'FILTER=\"PASS\"' -c 1 -O b -m2 -M2 --force-samples --types snps {input.bcf} -o {output.bcf}"

rule setid:
    input:
        vcf=rules.extract.output.bcf,
    output:
        #bcf=protected("tmp/{s_studyid}.{chrom}.annotated.bcf")
        bcf="tmp/{s_studyid}.{chrom}.annotated.bcf"
    conda: "../envs/bcftools.yaml"
    shell: "bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -O b -o {output} {input.vcf}"

rule convert_to_plink:
    input: bcf=rules.setid.output.bcf
    output:
        #bed=protected(TMP/"{s_studyid}_annotated_plink_chr{chrom}.bed"),
        #bim=protected(TMP/"{s_studyid}_annotated_plink_chr{chrom}.bim"),
        #fam=protected(TMP/"{s_studyid}_annotated_plink_chr{chrom}.fam"),
        bed=TMP/"{s_studyid}_annotated_plink_chr{chrom}.bed",
        bim=TMP/"{s_studyid}_annotated_plink_chr{chrom}.bim",
        fam=TMP/"{s_studyid}_annotated_plink_chr{chrom}.fam",
    conda: "../envs/plink2a.yaml"
    params: out=lambda w: TMP/f"{w.s_studyid}_annotated_plink_chr{w.chrom}"
    shell:
        """plink --allow-extra-chr  --bcf {input} --make-bed --out {params.out}"""
        
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
        bed=rules.annotate_PAR.output.bed,
        bim=rules.annotate_PAR.output.bed,
        fam=rules.annotate_PAR.output.bed,
    output:
        bed=TMP/"{s_studyid}_annotated_plink_merged.vcf.gz",
    conda: "../envs/plink2a.yaml"
    shell: """plink --bfile tmp/{wildcards.s_studyid}_annotated_plink_merged --out tmp/{wildcards.s_studyid}_annotated_plink_merged --recode vcf bgz"""

rule create_pedigree_fam_file:
    input: fam=rules.merge_and_filter.output.fam
    output: TMP/"{s_studyid}_annotated_plink_merged_ped.fam"
    conda: "../envs/r.yaml"
    script: "../scripts/R/put_in_pedigree.R"
    
rule king_duplicate_check:
    input:
        bed=TMP/"{s_studyid}_annotated_plink_merged_refformat.bed",
        bim=TMP/"{s_studyid}_annotated_plink_merged_refformat.bim",
        fam=TMP/"{s_studyid}_annotated_plink_merged_refformat.fam",
    output: temp(TMP/"{s_studyid}_king_duplicate.con")
    conda: "../envs/king.yaml"
    shell: "king -b {input.bed} --duplicate --prefix tmp/{wildcards.s_studyid}_king_duplicate"

rule collate_removed_samples:
    input:
        manifest="multiomics/{s_studyid}/ANNOTATED_MANIFEST.csv",
        flags=FLAGS/'{s_studyid}_samples.csv',
    output:
        txt="tmp/{s_studyid}/remove_samples.txt",
        csv="output_freezes/{s_studyid}/samples_removed.csv",
    script: "../scripts/python/collate_removed_samples.py"

rule filter:
    input:
        bed="tmp/{s_studyid}_annotated_plink_chr{chrom}.bed",
        bim="tmp/{s_studyid}_annotated_plink_chr{chrom}.bim",
        fam="tmp/{s_studyid}_annotated_plink_chr{chrom}.fam",
        rem="tmp/{s_studyid}/remove_samples.txt",
        exclude="flags/{s_studyid}_markers.csv",
    output:
        vcf=TMP/"{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.vcf.gz",
    params:
        out="tmp/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}"
    conda: "../envs/plink2a.yaml"
    shell: "plink --bed {input.bed} \
                  --bim {input.bim} \
                  --fam {input.fam} \
                  --remove {input.rem} \
                  --exclude {input.exclude} \
                  --out {params.out} \
                  --recode vcf bgz"

rule convert_and_publish:
    input: vcf=TMP/"{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.vcf.gz",
    output:
        bcf="output_freezes/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.bcf",
        csi="output_freezes/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.bcf.csi",
    conda: "../envs/bcftools.yaml"
    shell: """bcftools convert -O b -o {output.bcf} {input.vcf}; \
            bcftools index -f -c {output.bcf}"""

rule pep:
    input:
        manifest='multiomics/{s_studyid}/ANNOTATED_MANIFEST.csv',
        bcfs=expand('output_freezes/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.bcf', s_studyid=STUDIES, chrom=CHROMOSOMES),
        indexes=expand('output_freezes/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.bcf.csi', s_studyid=STUDIES, chrom=CHROMOSOMES),
    output:
        pep='output_freezes/{s_studyid}/{s_studyid}_freeze.10.pep.yaml',
        sample_table='output_freezes/{s_studyid}/{s_studyid}_freeze.10.pep.csv'
    conda: "../envs/pep.yaml"
    params:
        chromosomes=CHROMOSOMES,
        cmdline=CMDLINE,
        timestamp=TIMESTAMP,
        reference_genome='GRCh38'
    script: "../scripts/python/write_pep_files.py"

rule md5sums:
    input:
        bcfs=expand('output_freezes/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.bcf', s_studyid=STUDIES, chrom=CHROMOSOMES),
    output:
        md5sums='output_freezes/{s_studyid}/md5sums.txt',
    script: "md5sum {input.bcfs} > {output.md5sums}"

rule publish_all:
    input:
        bcfs=expand('output_freezes/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.bcf', s_studyid=STUDIES, chrom=CHROMOSOMES),
        indexes=expand('output_freezes/{s_studyid}/{s_studyid}_freeze.10_chr{chrom}.bcf.csi', s_studyid=STUDIES, chrom=CHROMOSOMES),


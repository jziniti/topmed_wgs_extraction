from snakemake.utils import min_version
min_version("6.0")

# TARGETS = expand(TMP/'GECOPD_annotated_plink_merged.{suffix}', suffix=('frq', 'imiss', 'lmiss', 'het', 'hwe', 'sexcheck'))
rule: input: TARGETS

module plink_commands:
    snakefile: "../cdnm/rules/plink.smk"

use rule * from plink_commands

""" use rule plink_freq from plink_commands with:
    input:
        bed='tmp/{s_studyid}_annotated_plink_merged.bed',
        bim='tmp/{s_studyid}_annotated_plink_merged.bim',
        fam='tmp/{s_studyid}_annotated_plink_merged.fam',

use rule plink_missing from plink_commands with:
    input:
        bed='tmp/{s_studyid}_annotated_plink_merged.bed',
        bim='tmp/{s_studyid}_annotated_plink_merged.bim',
        fam='tmp/{s_studyid}_annotated_plink_merged.fam',

use rule plink_het from plink_commands with:
    input:
        bed='tmp/{s_studyid}_annotated_plink_merged.bed',
        bim='tmp/{s_studyid}_annotated_plink_merged.bim',
        fam='tmp/{s_studyid}_annotated_plink_merged.fam',

use rule plink_hwe from plink_commands with:
    input:
        bed='tmp/{s_studyid}_annotated_plink_merged.bed',
        bim='tmp/{s_studyid}_annotated_plink_merged.bim',
        fam='tmp/{s_studyid}_annotated_plink_merged.fam',

use rule plink_sexcheck from plink_commands with:
    input:
        bed='tmp/{s_studyid}_annotated_plink_merged.bed',
        bim='tmp/{s_studyid}_annotated_plink_merged.bim',
        fam='tmp/{s_studyid}_annotated_plink_merged.fam',
 """
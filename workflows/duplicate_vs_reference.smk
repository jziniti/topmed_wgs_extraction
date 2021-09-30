from snakemake.utils import min_version
min_version("6.0")

TMP = Path('tmp')

rule convert_to_chrpos:
    input:
        bed=TMP/"{s_studyid}_annotated_plink_merged.bed",
        bim=TMP/"{s_studyid}_convert_chrpos.bim",
        fam=TMP/"{s_studyid}_annotated_plink_merged.fam",
    output:
        bed=TMP/"{s_studyid}_annotated_plink_merged_chrpos.bed",
        bim=TMP/"{s_studyid}_annotated_plink_merged_chrpos.bim",
        fam=TMP/"{s_studyid}_annotated_plink_merged_chrpos.fam",
    params: prefix=lambda w: TMP/f"{w.s_studyid}_annotated_plink_merged_chrpos",
    conda: "../envs/plink2a.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out {params.prefix} --make-bed"

module king:
    snakefile: "../king/rules/duplicate_vs_reference.smk"

use rule duplicate_vs_reference from king with:
    input:
        qbed=TMP/"{s_studyid}_annotated_plink_merged_chrpos.bed",
        qbim=TMP/"{s_studyid}_annotated_plink_merged_chrpos.bim",
        qfam=TMP/"{s_studyid}_annotated_plink_merged_chrpos.fam",
        rbed=lambda w: f"{config[w['s_studyid']]['known_good_reference']}.bed",
        rbim=lambda w: f"{config[w['s_studyid']]['known_good_reference']}.bim",
        rfam=lambda w: f"{config[w['s_studyid']]['known_good_reference']}.fam",
    output:
        kin=TMP/"{s_studyid}_reference_concordance.con",
    params:
        Q=lambda w: TMP/f"{w.s_studyid}_annotated_plink_merged_chrpos",
        R=lambda w: f"{config[w['s_studyid']]['known_good_reference']}",
        prefix=lambda w: TMP/f"{w.s_studyid}_reference_concordance",

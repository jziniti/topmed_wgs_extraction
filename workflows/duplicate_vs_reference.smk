from snakemake.utils import min_version
min_version("6.0")

TMP = Path('tmp')

rule convert_bim_to_civic:
    input: bim=TMP/"{s_studyid}_annotated_plink_merged.bim",
    output: bim=temp(TMP/"{s_studyid}_annotated_plink_merged_civic.bim")
    shell: "awk '{{print $1,$1\":\"$4\":\"$5\":\"$6,$3,$4,$5,$6}}' {input.bim} > {output.bim}"


rule convert_bim_to_chrpos:
    input: bim=TMP/"{s_studyid}_annotated_plink_merged.bim"
    output: bim=temp(TMP/"{s_studyid}_annotated_plink_merged_chrpos.bim")
    shell: "awk '{{print $1,$1\":\"$4,$3,$4,$5,$6}}' {input.bim} > {output.bim}"


rule match_format_to_reference:
    input:
        bed=TMP/"{s_studyid}_annotated_plink_merged.bed",
        bim=lambda w: TMP/f"{w.s_studyid}_annotated_plink_merged_{config[w.s_studyid].get('refformat', 'civic')}.bim",
        fam=TMP/"{s_studyid}_annotated_plink_merged.fam",
    output:
        bed=TMP/"{s_studyid}_annotated_plink_merged_refformat.bed",
        bim=TMP/"{s_studyid}_annotated_plink_merged_refformat.bim",
        fam=TMP/"{s_studyid}_annotated_plink_merged_refformat.fam",
    params: prefix=lambda w: TMP/f"{w.s_studyid}_annotated_plink_merged_refformat",  ### This *must* match what is in output
    conda: "../envs/plink2a.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out {params.prefix} --make-bed"

module king:
    snakefile: "../king/rules/duplicate_vs_reference.smk"

### FIXME ... this rule is hard to generalize ... the reference datasets
use rule duplicate_vs_reference from king with:
    input:
        qbed=rules.match_format_to_reference.output.bed,
        qbim=rules.match_format_to_reference.output.bim,
        qfam=rules.match_format_to_reference.output.fam,
        rbed=lambda w: f"{config[w['s_studyid']]['known_good_reference']}.bed",
        rbim=lambda w: f"{config[w['s_studyid']]['known_good_reference']}.bim",
        rfam=lambda w: f"{config[w['s_studyid']]['known_good_reference']}.fam",
    output:
        kin=TMP/"{s_studyid}_reference_concordance.con",
    resources:
        mem_free="200G",
        mem_mb="2000",
    params:
        Q=lambda w: TMP/f"{w.s_studyid}_annotated_plink_merged_chrpos",
        R=lambda w: f"{config[w['s_studyid']]['known_good_reference']}",
        prefix=lambda w: TMP/f"{w.s_studyid}_reference_concordance",

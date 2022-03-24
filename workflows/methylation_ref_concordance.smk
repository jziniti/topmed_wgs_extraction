from snakemake.utils import min_version
min_version("6.0")

TMP = Path('tmp')
OUTPUT_BASE = 'methylation_ref_concordance/gecopd_freeze.10_concordance_report'

rule: input: [TMP/'GECOPD_methylation_unified.kin', TMP/'GECOPD_methylation_freeze10.kin']

# N_BATCHES = config['GECOPD']['batch_count']

rule convert_chrpos:
    input:
        bim=TMP/"{s_studyid}_annotated_plink_merged.bim"
    output:
        bim=temp(TMP/"{s_studyid}_annotated_convert_chrpos.bim")
    shell: "awk '{{print $1,$1\":\"$4,$3,$4,$5,$6}}' {input.bim} > {output.bim}"

rule convert_to_pedmap:
    input:
        bed=TMP/"{s_studyid}_annotated_plink_merged.bed",
        bim=TMP/"{s_studyid}_annotated_convert_chrpos.bim",
        fam=TMP/"{s_studyid}_annotated_plink_merged.fam",
        extract=TMP/"{s_studyid}_extract.txt"
    output:
        ped=TMP/"{s_studyid}_annotated_plink_merged.ped",
        map=TMP/"{s_studyid}_annotated_plink_merged.map",
    params: out="tmp/{s_studyid}_annotated_plink_merged",
    conda: "../envs/plink.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out {params.out} --extract {input.extract} --recode"

module methylation_ref_concordance:
    snakefile: "../850K_qc_pipeline/workflows/external-reference-concordance.smk"
    config: config[list(config.keys())[0]]
    # config: config

use rule * from methylation_ref_concordance as mrc_* 

use rule write_extract_file from methylation_ref_concordance as mrc_write_extract_file with:
    output: temp(TMP/'{s_studyid}_extract.txt')

#use rule compare_genotypes from methylation_ref_concordance as mrc_compare_genotypes with:
#    input:
#        control_snp_calls=expand("/proj/regeps/regep00/studies/GECOPD/data/methylation/TopMed/data/freezes/20200706/LEVEL2/control_snp_genotypes.{idx}.txt", idx=range(N_BATCHES)),
#        ped=f'{OUTPUT_BASE}.ped',
#        map=f'{OUTPUT_BASE}.map',
#        assumed=srcdir("../850K_qc_pipeline/assume.txt"),
#        snpmap=srcdir("../850K_qc_pipeline/epic_control_snp_list_GRCh38.p7.txt"),
#    output:
#        kin=temp(TMP/'GECOPD_methylation_unified.kin'),
#        kin0=temp(TMP/'GECOPD_methylation_unified.kin0'),

use rule compare_genotypes from methylation_ref_concordance as mfc_compare_genotypes with:
    input:
        control_snp_calls=lambda w: expand(Path(config[w.s_studyid]['methylation_freeze_path'])/f"control_snp_genotypes.{{idx}}.txt", idx=range(config[w.s_studyid]['methylation_batch_count'])),
        ped='tmp/{s_studyid}_annotated_plink_merged.ped',
        map='tmp/{s_studyid}_annotated_plink_merged.map',
        assumed=srcdir("../850K_qc_pipeline/assume.txt"),
        snpmap=srcdir("../850K_qc_pipeline/epic_control_snp_list_GRCh38.p7.txt"),
    output:
        kin=temp(TMP/'{s_studyid}_methylation_freeze10.kin'),
        kin0=temp(TMP/'{s_studyid}_methylation_freeze10.kin0'),
from snakemake.utils import min_version
min_version("6.0")

TMP = Path('tmp')
OUTPUT_BASE = 'methylation_ref_concordance/gecopd_freeze.10_concordance_report'

rule: input: [TMP/'GECOPD_methylation_unified.kin', TMP/'GECOPD_methylation_freeze10.kin']

N_BATCHES = config['GECOPD']['batch_count']

rule convert_to_pedmap:
    input:
        bed="tmp/GECOPD_annotated_plink_merged.bed",
        bim="tmp/GECOPD_annotated_plink_merged.bim",
        fam="tmp/GECOPD_annotated_plink_merged.fam",
        extract="tmp/extract.txt"
    output:
        ped="tmp/GECOPD_annotated_plink_merged.ped",
        map="tmp/GECOPD_annotated_plink_merged.map",
    params: out="tmp/GECOPD_annotated_plink_merged",
    conda: "../envs/plink.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out {params.out} --extract {input.extract} --recode"

module methylation_ref_concordance:
    snakefile: "../850K_qc_pipeline/workflows/external-reference-concordance.smk"
    config: config['GECOPD']

use rule * from methylation_ref_concordance as mrc_* 

use rule compare_genotypes from methylation_ref_concordance as mrc_compare_genotypes with:
    input:
        control_snp_calls=expand("/proj/regeps/regep00/studies/GECOPD/data/methylation/TopMed/data/freezes/20200706/LEVEL2/control_snp_genotypes.{idx}.txt", idx=range(N_BATCHES)),
        ped=f'{OUTPUT_BASE}.ped',
        map=f'{OUTPUT_BASE}.map',
        assumed=srcdir("../850K_qc_pipeline/assume.txt"),
        snpmap=srcdir("../850K_qc_pipeline/epic_control_snp_list_GRCh38.p7.txt"),
    output:
        kin=temp(TMP/'GECOPD_methylation_unified.kin'),
        kin0=temp(TMP/'GECOPD_methylation_unified.kin0'),

use rule compare_genotypes from methylation_ref_concordance as mfc_compare_genotypes with:
    input:
        control_snp_calls=expand("/proj/regeps/regep00/studies/GECOPD/data/methylation/TopMed/data/freezes/20200706/LEVEL2/control_snp_genotypes.{idx}.txt", idx=range(N_BATCHES)),
        ped=f'tmp/GECOPD_annotated_plink_merged.ped',
        map=f'tmp/GECOPD_annotated_plink_merged.map',
        assumed=srcdir("../850K_qc_pipeline/assume.txt"),
        snpmap=srcdir("../850K_qc_pipeline/epic_control_snp_list_GRCh38.p7.txt"),
    output:
        kin=temp(TMP/'GECOPD_methylation_freeze10.kin'),
        kin0=temp(TMP/'GECOPD_methylation_freeze10.kin0'),
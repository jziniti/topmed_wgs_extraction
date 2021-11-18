
TARGETS.append('plate103/ANNOTATED_MANIFEST.csv')
TARGETS.append('plate103/filter.1/ANNOTATED_MANIFEST.csv')

rule: input: TARGETS

rule p103_filter_annotated_manifest:
    input:
        manifest="multiomics/GECOPD/ANNOTATED_MANIFEST.csv",
        plate_manifests=['plate103/PL-01007365_v2.csv', 'plate103/PL-01007396_v2.csv', 'plate103/PL-01007962_v2.csv'],
        sample_metadata="plate103/p103_sapphire8_sample_metadata.csv"
    output: manifest="plate103/ANNOTATED_MANIFEST.csv"
    script: "../scripts/python/plate103_filter_manifest.py"

rule p103_1_filter_annotated_manifest:
    input:
        manifest="multiomics/GECOPD.1/ANNOTATED_MANIFEST.csv",
        plate_manifests=['plate103/PL-01007365_v2.csv', 'plate103/PL-01007396_v2.csv', 'plate103/PL-01007962_v2.csv'],
        sample_metadata="plate103/p103_sapphire8_sample_metadata.csv"
    output: manifest="plate103/filter.1/ANNOTATED_MANIFEST.csv"
    script: "../scripts/python/plate103_filter_manifest.py"

rule p103_sapphire8_sample_metadata:
    input: plate_manifests=['plate103/PL-01007365_v2.csv', 'plate103/PL-01007396_v2.csv', 'plate103/PL-01007962_v2.csv'],
    output: "plate103/p103_sapphire8_sample_metadata.csv"
    conda: "../cdnm/envs/pandas_oracle.yaml"
    script: "../scripts/python/get_sapphire8_sample_metadata.py"

### FIXME: write rule for plate103 drop-list"
rule p103_make_drop_list_1:
    input:
        manifest="plate103/ANNOTATED_MANIFEST.csv",
        manual_drop_list="plate103/drop.manual.txt",
        plate103_pre_drop_list="plate103/NWDID_DROP_Plate_103.txt"
    output: "plate103/drop.txt"
    shell: "cat {input.manual_drop_list} {input.plate103_pre_drop_list}| sort | uniq | awk '{{print \"0\\t\"$1}}' > {output}"

rule p103_remove_bad_samples:
    input:
        bed=TMP/"GECOPD_annotated_plink_merged.bed",
        bim=TMP/"GECOPD_annotated_plink_merged.bim",
        fam=TMP/"GECOPD_annotated_plink_merged.fam",
        droplist="plate103/drop.txt"
    output:
        bed=TMP/"GECOPD_annotated_plink_merged_1.bed",
        bim=TMP/"GECOPD_annotated_plink_merged_1.bim",
        fam=TMP/"GECOPD_annotated_plink_merged_1.fam",
    params:
        prefix=TMP/"GECOPD_annotated_plink_merged_1"
    conda: "../cdnm/envs/plink.yaml"
    shell:  "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out {params.prefix} --remove {input.droplist} --make-bed"

#module multiomics: 
#    snakefile: "../topmed_multiomics_qc/workflows/multiomics_qc.smk"
#    config: "conf/plate103.filter1.yaml"

#use rule * from multiomics as p103_*

use rule add_manual_flags from multiomics as p103_add_manual_flags with:
    input:
        manifest=lambda w: TMP/f'GECOPD.1/MANIFEST.flag_sample_groups.csv'
    output:
        df='multiomics/GECOPD.1/ANNOTATED_MANIFEST.csv',



TARGETS.append('plate103/ANNOTATED_MANIFEST.csv')

rule: input: TARGETS

rule p103_filter_annotated_manifest:
    input:
        manifest="multiomics/GECOPD/ANNOTATED_MANIFEST.csv",
        plate_manifests=['plate103/PL-01007365_v2.csv', 'plate103/PL-01007396_v2.csv', 'plate103/PL-01007962_v2.csv'],
        sample_metadata="plate103/p103_sapphire8_sample_metadata.csv"
    output: manifest="plate103_ANNOTATED_MANIFEST.csv"
    script: "../scripts/python/plate103_filter_manifest.py"

rule p103_sapphire8_sample_metadata:
    input: plate_manifests=['plate103/PL-01007365_v2.csv', 'plate103/PL-01007396_v2.csv', 'plate103/PL-01007962_v2.csv'],
    output: "plate103/p103_sapphire8_sample_metadata.csv"
    conda: "../cdnm/envs/pandas_oracle.yaml"
    script: "../scripts/python/get_sapphire8_sample_metadata.py"
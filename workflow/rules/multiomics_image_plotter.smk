MANIFEST_PATH = config.get("manifest_path", "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm/plate103/ANNOTATED_MANIFEST.csv")
C_DATA_PATH = config.get("c_data_path", "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm/plate103/c_data.build_sample_groups")
GROUP_IDS = config.get("draw_concordance_groups", ["all"])
SHAPE_BY = config.get("shape_by", "ObservedGender")

rule: input: expand('images/group{group_id}.png', group_id=GROUP_IDS)

rule render_concordance_group:
    input:
        annotated_manifest=MANIFEST_PATH,
        concordance_data=C_DATA_PATH
    output:
        expand('images/group{group_id}.png', group_id=GROUP_IDS)
    params:
        group_ids=GROUP_IDS,
        shape_by=SHAPE_BY
    conda: "../envs/cdnm-jupyter-python-3.7.6.yaml"
    script: "../scripts/python/image_plotter.py"


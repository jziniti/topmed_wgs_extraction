
GROUP_IDS = config.get('draw_concordance_groups', [])
C_DATA_PATH = config.get("c_data_path", "plate103/c_data.build_sample_groups")
SAMPLE_GROUP_PATH = config.get("sample_group_path", "plate103/sample_groups.flag_sample_groups")
MANIFEST_PATH = config.get("manifest_path", "plate103/ANNOTATED_MANIFEST.csv")

TARGETS = expand('images/concordance_groups/{group_id}.png', group_id=GROUP_IDS)

rule: input: TARGETS

rule render_concordance_group:
    input:
        concordance_data=C_DATA_PATH,
        sample_group_metadata=SAMPLE_GROUP_PATH,
        annotated_manifest=MANIFEST_PATH,
    output:
        expand('images/concordance_groups/{group_id}.png', group_id=GROUP_IDS)
    params:
        group_ids=GROUP_IDS
    conda: "../cdnm/envs/r-visnetwork.yaml"
    script: "../scripts/R/render_concordance_group.R"

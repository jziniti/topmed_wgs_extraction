
TARGETS = expand('images/concordance_groups/{group_id}.png', group_id=config.get('draw_concordance_groups', []))

rule: input: TARGETS

rule render_concordance_group:
    input:
        concordance_data="plate103/c_data.build_sample_groups",
        sample_group_metadata="plate103/sample_groups.flag_sample_groups",
        annotated_manifest="plate103/ANNOTATED_MANIFEST.csv"
    output:
        expand('images/concordance_groups/{group_id}.png', group_id=config.get('draw_concordance_groups', []))
    params:
        group_ids=config.get('draw_concordance_groups', [])
    conda: "../cdnm/envs/r-visnetwork.yaml"
    script: "../scripts/R/render_concordance_group.R"
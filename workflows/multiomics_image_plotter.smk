SAMPLE_GROUP_PATH="/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm/plate103/ANNOTATED_MANIFEST.csv"
my_group_id = 10

rule: input: expand('images/{group_id}.txt', group_id=my_group_id)

rule render_concordance_group:
    input:
        concordance_data=SAMPLE_GROUP_PATH,
    output:
        expand('images/{group_id}.txt', group_id=my_group_id)
    params:
        group_ids=my_group_id
    conda: "../envs/cdnm-jupyter-python-3.7.6.yaml"
    script: "../scripts/python/image_plotter.py"
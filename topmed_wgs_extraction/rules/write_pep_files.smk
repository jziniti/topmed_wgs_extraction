rule write_pep_files:
    input:
        sample_manifest='multiomics/{configkey}/ANNOTATED_MANIFEST.csv',
    output:
        project_metadata='{configkey}_pep.yaml',
        sample_table='{configkey}_samples_pep.csv',
    conda: "../envs/pep.yaml"
    params: config=config
    script: "../scripts/python/write_pep_files.py"
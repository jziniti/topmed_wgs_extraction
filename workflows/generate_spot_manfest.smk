TARGETS = expand('spot/{s_studyid}/MANIFEST.csv', s_studyid=(config.keys()))

rule: input: TARGETS

rule generate_spot_manifest:
    input: manifest="multiomics/{s_studyid}/ANNOTATED_MANIFEST.csv"
    output: manifest="spot/{s_studyid}/MANIFEST.csv"
    script: "../scripts/python/generate_spot_manifest.py"
    

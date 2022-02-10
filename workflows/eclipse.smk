rule hrc_concordance:
    input:
        qbed=TMP/"ECLPSE_annotated_plink_merged.bed",
        rbed=lambda w: f"{config['ECLPSE']['known_good_reference']}.bed",
    output:
        csv=TMP/"ECLPSE_hrcimpute_concordance.csv"
    params:
        qbed=TMP/"ECLPSE_annotated_plink_merged",
        rbed=config['ECLPSE']['known_good_reference'],
    conda: "../envs/full-similarity-matrix.yaml"
    shell: "cdnm-wf king-similarity-matrix \
                        --query-bed {params.qbed}\
                        --reference-bed {params.rbed}\
                        --output {output.csv}"
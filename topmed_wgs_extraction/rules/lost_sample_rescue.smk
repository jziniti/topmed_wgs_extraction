rule run_lost_sample_qc:
    input:
    output: "tmp/lost_samples.out"
    conda: "../envs/r.yaml"
    script: "../notebooks/lost_samples.Rmd"

rule king:


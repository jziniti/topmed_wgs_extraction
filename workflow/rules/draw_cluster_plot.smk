rule draw_cluster_plot:
    input:
        manifest="multiomics/GECOPD/ANNOTATED_MANIFEST.csv",
        c_data="multiomics/GECOPD/c_data.build_sample_groups",
    output:
    conda: "../envs/pydot.yaml"
    script: "../scripts/python/draw_cluster_plot.py"
            

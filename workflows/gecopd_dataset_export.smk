rule gecopd_collate_marker_flags:
    input:
        flags="flags/GECOPD_markers.csv"
    output:
        flags="flags/GECOPD_markers.tsv"

rule reassign_samples:
     input:
        plate103="plate103/reassign_samples.csv",
     output:
     shell:

rule gecopd_collate_sample_flags:
    input:
        flags="flags/GECOPD_markers.csv",
        plate103="plate103/drop_samples.csv",
    output:
        flags="flags/GECOPD_markers.tsv"
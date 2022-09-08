from pathlib import Path

TMP = Path('tmp')
TARGETS = []
INPUT_PATHS = {}

with open('conf/all_rna_batches.txt') as fh:
    for line in fh:
        batch = line.strip().split('/')[-3]
        TARGETS.append(TMP/f"{batch}_king.csv")
        INPUT_PATHS[batch] = line.strip()

rule: input: TARGETS 

rule generate_config:
    input:
        qbed=lambda w: INPUT_PATHS[w.batch],
        rbed=TMP/"GECOPD_annotated_plink_merged.bed",
    output: TMP/"{batch}.yaml"
    params:
        qbase=lambda w: INPUT_PATHS[w.batch].replace(".bed", ""),
        rbase=TMP/"GECOPD_annotated_plink_merged",
    run:
        with open(output[0], 'w') as fh:
            fh.write(f'query_bed: "{input.qbed}"\n')
            fh.write(f'reference_bed: "{input.rbed}"\n')
            fh.write(f'prefix: "tmp/{wildcards.batch}_king"\n')

rule run_king_similarity_matrix:
    input:
        qbed=lambda w: INPUT_PATHS[w.batch],
        rbed=TMP/"GECOPD_annotated_plink_merged.bed",
        config=TMP/"{batch}.yaml",
    output: TMP/"{batch}_king.csv"
    conda: "../envs/full-similarity-matrix.yaml"
    shell: "cdnm-wf king-similarity-matrix --config {input.config}"

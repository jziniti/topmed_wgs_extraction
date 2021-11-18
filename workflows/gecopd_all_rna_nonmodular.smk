from pathlib import Path

TMP = Path('tmp')
TARGETS = []
INPUT_PATHS = {}

with open('conf/all_rna_batches.txt') as fh:
    for line in fh:
        batch = line.strip().split('/')[-3]
        TARGETS.append(TMP/f"{batch}.con")
        INPUT_PATHS[batch] = line.strip()

wildcard_constraints:
    batch='batch_.*'

rule: input: TARGETS 

rule convert_to_chrpos:
    input: lambda w: f'{INPUT_PATHS[w.batch].replace("bed", "bim")}'
    output: TMP/'{batch}_chrpos.bim'
    shell: "cat {input} | awk '{{ print $1,$1\":\"$4,$3,$4,$5,$6}}' > {output}"

rule extract_keep_list:
    input: lambda w: f'{INPUT_PATHS[w.batch].replace("bed", "fam")}'
    output: TMP/'{batch}.keep'
    shell: "grep 'S-' {input} > {output}"

rule prepare_for_king:
    input:
        bed=lambda w: f'{INPUT_PATHS[w.batch]}',
        bim=TMP/'{batch}_chrpos.bim',
        keep=TMP/'{batch}.keep',
        fam=lambda w: f'{INPUT_PATHS[w.batch].replace("bed", "fam")}',
    output: TMP/"{batch}.bed"
    conda: "../envs/plink.yaml"
    params: out=lambda w: TMP/f"{w.batch}",
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --keep {input.keep} --out {params.out} --make-bed"

rule run_king_similarity_matrix:
    input:
        qbed=TMP/"{batch}.bed",
        rbed=TMP/"GECOPD_annotated_plink_merged_chrpos.bed",
        config=TMP/"{batch}_config.yaml",
    output: TMP/"{batch}.con"
    params: prefix=lambda w: TMP/f'{w.batch}'
    conda: "../envs/king.yaml"
    shell: "king -b {input.qbed},{input.rbed} --prefix {params.prefix}"

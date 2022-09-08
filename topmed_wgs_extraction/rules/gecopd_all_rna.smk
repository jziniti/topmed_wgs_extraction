from pathlib import Path

TMP = Path('tmp')
TARGETS = [TMP/"GECOPD_all_rna_manifest.txt"]

INPUT_PATHS = {}
BATCH_OUTPUTS = []

with open('conf/all_rna_batches.txt') as fh:
    for line in fh:
        batch = line.strip().split('/')[-3]
        BATCH_OUTPUTS.append(TMP/f"{batch}_king.csv")
        INPUT_PATHS[batch] = line.strip()

rule: input: TARGETS 

wildcard_constraints:
    batch='batch_.*'

rule extract_keep_list:
    input: lambda w: f'{INPUT_PATHS[w.batch].replace("bed", "fam")}'
    output: TMP/'{batch}.keep'
    shell: "grep 'S-' {input} > {output}"

rule prepare_for_king:
    input:
        bed=lambda w: f'{INPUT_PATHS[w.batch]}',
        bim=lambda w: f'{INPUT_PATHS[w.batch].replace("bed", "bim")}',
        fam=lambda w: f'{INPUT_PATHS[w.batch].replace("bed", "fam")}',
        keep=TMP/'{batch}.keep',
    output: TMP/"{batch}.bed"
    conda: "../envs/plink.yaml"
    params: out=lambda w: TMP/f"{w.batch}",
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --keep {input.keep} --out {params.out} --make-bed"

rule generate_config:
    input:
        qbed=TMP/"{batch}.bed",
        rbed=TMP/"GECOPD_annotated_plink_merged.bed",
    output: TMP/"{batch}.yaml"
    params:
        rbase=TMP/"GECOPD_annotated_plink_merged",
    run:
        with open(output[0], 'w') as fh:
            fh.write(f'query_bed: "tmp/{wildcards.batch}.bed"\n')
            fh.write(f'reference_bed: "{params.rbase}.bed"\n')
            fh.write(f'prefix: "tmp/{wildcards.batch}_king.csv"\n')
            #fh.write(f'dataset_base: "tmp/{wildcards.batch}"\n')
            #fh.write(f'reference_base: "{params.rbase}"\n')
            #fh.write(f'output_file: "tmp/{wildcards.batch}_king.csv"\n')

rule run_king_similarity_matrix:
    input:
        qbed=TMP/"{batch}.bed",
        rbed=TMP/"GECOPD_annotated_plink_merged.bed",
        config=TMP/"{batch}.yaml",
    output: TMP/"{batch}_king.csv"
    conda: "../envs/full-similarity-matrix.yaml"
    shell: "cdnm-wf king-similarity-matrix --config {input.config}"

rule collate_rna_wgs_batch_concordance:
    input:
        batch_king_results=expand(TMP/"{batch}_king_rna_wgs_filtered.csv", batch=INPUT_PATHS.keys()),
    output:
        rna_wgs_king_results=TMP/"GECOPD_rna_wgs_king_results.csv",
    shell: "(head -n 1 {input.batch_king_results[0]}; grep --no-filename \"S-\" {input.batch_king_results}) > {output}"
    #shell: "(head -n 1 {input.batch_king_results[0]}; tail -n +2 {input.batch_king_results}) > {output}"

#rule collate
#    input:
#        batch_king_results=BATCH_OUTPUTS,
#    output:
#        TMP/"GECOPD_rna_rna_collated_king_results.csv",
#    shell:
#        ""

rule adjust_concordance_data_to_masterfile:
    input:
        masterfile='/proj/edith/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/archivalroot/Masterfile/master.file.freeze4.txt',
        batch=TMP/"{batch}_king.csv",
    output:
        csv=TMP/"{batch}_king_rna_wgs_filtered.csv",
    script: "../scripts/python/filter_king_results_rna_wgs.py"

rule read_rna_masterfile:
    input:
        masterfile='/proj/edith/regeps/regep00/studies/COPDGene/data/rnaseq/BLOOD/archivalroot/Masterfile/master.file.freeze4.txt',
        ## phenofile='/proj/regeps/regep00/studies/COPDGene/data/pheno/Final10000_Dataset_12MAR13.txt',
        phenofile='/proj/regeps/regep00/studies/COPDGene/data/pheno/COPDGene_P1P2P3_Flat_SM_NS_Nov21.txt',
    output:
        manifest=TMP/"GECOPD_all_rna_manifest.txt",
        gobs="metadata/gecopd_rna_fake_gender.csv",
        reassignments="metadata/GECOPD_rna_reassign.csv",
    script: "../scripts/python/filter_rna_masterfile.py"

rule merge_resta_king_results:
    input:
        manifest="tmp/GECOPD_all_rna_manifest.txt",
        king_results = [
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_201221_NB551475_0238_AHNKL7BGXG_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_201223_NB552221_0052_AH3FJFBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_201228_NB552221_0053_AH3G3VBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210114_NB552221_0059_AHC3HLBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210129_NB552221_0061_AHC53FBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210202_NB551475_0245_AHGNFWBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210202_NB552221_0062_AHGNH2BGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210204_NB551475_0246_AHC5H7BGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210204_NB552221_0063_AHH3J5BGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210206_NB552221_0064_AHHGC5BGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210208_NB552221_0065_AHGN2TBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210209_NB551475_0249_AHGVGWBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210301_NB551475_0252_AHNCT3BGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210301_NB552221_0069_AHNFFMBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210302_NB551475_0253_AHNCVYBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210302_NB552221_0070_AHNFFTBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210318_NB551475_0256_AHNFFHBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_210318_NB552221_0073_AHNG7GBGXH_grch38_gtf/king.csv',
        '/udd/resta/COPDGene_rna_mrna_star_alignment_qc/rnaseq_qc_app/data/COPDGene_211007_NB552221_0097_AH53HHBGXK_grch38_gtf/king.csv',
        ]
    output: df="tmp/GECOPD_rna_wgs_resta_king_results.csv"
    script: "../scripts/python/merge_resta_king_results.py"


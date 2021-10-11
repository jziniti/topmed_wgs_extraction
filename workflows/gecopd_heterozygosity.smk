rule write_race_keep_file:
    input:
        fam=TMP/"GECOPD_annotated_plink_merged.fam",
        pheno="/proj/regeps/regep00/studies/COPDGene/data/pheno/Final10000_Phase1_Rev_28oct16.txt",
        stashq_file = lambda w: TMP/f"wgs.{SOURCE_FREEZE_FOR_STUDY.get('GECOPD', '10a')}.stashq.txt",
        idmap_file="/proj/regeps/regep00/studies/COPDGene/metadata/GECOPD_aliases.txt",
    output: keep=TMP/"GECOPD_{race_code}.keep"
    script: "../scripts/python/write_race_keep_file.py"

rule keep_by_race:
    input:
        bed=TMP/"GECOPD_annotated_plink_merged.bed",
        bim=TMP/"GECOPD_annotated_plink_merged.bim",
        fam=TMP/"GECOPD_annotated_plink_merged.fam",
        keep=TMP/"GECOPD_{race_code}.keep"
    output:
        bed=TMP/"GECOPD_{race_code,\d+}.bed",
        bim=TMP/"GECOPD_{race_code,\d+}.bim",
        fam=TMP/"GECOPD_{race_code,\d+}.fam",
    conda: "../cdnm/envs/plink.yaml"
    shell: "plink --bed {input.bed} -bim {input.bim} --fam {input.fam} --keep {input.keep} --out tmp/GECOPD_{wildcards.race_code} --make-bed"

rule merge_split_sexcheck_results:
    input:
        sexcheck_aa='tmp/GECOPD_2.sexcheck',
        sexcheck_nhw='tmp/GECOPD_1.sexcheck',
    output: merged_sexcheck='tmp/GECOPD_annotated_plink_merged.sexcheck'
    shell: "head -n 1 {input.sexcheck_aa} > {output}; grep --no-filename -v SNPSEX {input} >> {output}"

rule merge_split_het_results:
    input:
        aa='tmp/GECOPD_2.het',
        nhw='tmp/GECOPD_1.het',
    output: merged_sexcheck='tmp/GECOPD_annotated_plink_merged.het'
    shell: "head -n 1 {input.aa} > {output}; grep --no-filename -v FID {input} >> {output}"

rule hetereozygosity_notebook:
    input:
        pop_file = TMP/"{s_studyid}_annotated_plink_merged.fam",
        sexcheck_file= TMP/'{s_studyid}_annotated_plink_merged.sexcheck',
        hwe_files = expand(TMP/"GECOPD_{race_code}.hwe", race_code=[1,2]),
        het_files = expand(TMP/"GECOPD_{race_code}.het", race_code=[1,2]),
        stashq_file = lambda w: TMP/f"wgs.{SOURCE_FREEZE_FOR_STUDY.get(w.s_studyid, '10a')}.stashq.txt"
    output:
    log: notebook="notebook_logs/GECOPD_heterozygosity.ipynb"
    params:
    conda: "../envs/cdnm-jupyter-python-3.7.6.yaml"
    notebook: "../notebooks/standard_prelim_qc.ipynb"

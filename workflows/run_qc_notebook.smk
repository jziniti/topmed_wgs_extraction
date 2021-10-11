rule run_qc_notebook:
    input:
        pop_file = TMP/"{s_studyid}_annotated_plink_merged.fam",
        chrom_file = TMP/"{s_studyid}_annotated_plink_merged.bim",
        kin_file = TMP/'{s_studyid}_king_duplicate.con',
        sexcheck_file= TMP/'{s_studyid}_annotated_plink_merged.sexcheck',
        hwe_file = TMP/'{s_studyid}_annotated_plink_merged.hwe',
        frq_file = TMP/'{s_studyid}_annotated_plink_merged.frq',
        het_file = TMP/'{s_studyid}_annotated_plink_merged.het',
        imiss_file = TMP/'{s_studyid}_annotated_plink_merged.imiss',
        lmiss_file = TMP/'{s_studyid}_annotated_plink_merged.lmiss',
        stashq_file = lambda w: TMP/f"wgs.{SOURCE_FREEZE_FOR_STUDY.get(w.s_studyid, '10a')}.stashq.txt"
    output:
        sample_flags="flags/{s_studyid}_samples.csv",
        marker_flags="flags/{s_studyid}_markers.csv",
    log: notebook="notebook_logs/{s_studyid}_prelim_qc.ipynb"
    params:
        P=0.000001,
        f_bound=0.2,
        imiss_range=[10,20],
        lmiss_range=[1,5],
    conda: "../envs/cdnm-jupyter-python-3.7.6.yaml"
    notebook: "../notebooks/standard_prelim_qc.ipynb"

rule reaka_notebook:
    input:
        pop_file = TMP/"{s_studyid}_annotated_plink_merged.fam",
        chrom_file = TMP/"{s_studyid}_annotated_plink_merged.bim",
        kin_file = TMP/'tmp/{s_studyid}_king_duplicate.con',
        sexcheck_file= TMP/'{s_studyid}_annotated_plink_merged.sexcheck',
        hwe_file = TMP/'{s_studyid}_annotated_plink_merged.hwe',
        frq_file = TMP/'{s_studyid}_annotated_plink_merged.frq',
        het_file = TMP/'{s_studyid}_annotated_plink_merged.het',
        imiss_file = TMP/'{s_studyid}_annotated_plink_merged.imiss',
        lmiss_file = TMP/'{s_studyid}_annotated_plink_merged.lmiss',
        ## kin0_file = TMP/"{s_studyid}_king.kin0"
    output:
        sample_flags="flags/{s_studyid}_samples.csv",
        marker_flags="flags/{s_studyid}_markers.csv",
    log: notebook="notebook_logs/{s_studyid}_prelim_qc.ipynb"
    params:
        Error_one=1,
        Error_two=0.5,
        Kinship=0.354,
        SNPSEX=0,
        P=0.000001,
        f_bound=0.2
    conda: "../envs/cdnm-jupyter-python-3.7.6.yaml"
    notebook: "../notebooks/GECOPDNotebook1.ipynb"

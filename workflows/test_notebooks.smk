TARGETS = ['tmp/flags_remra.txt', 'tmp/flags_reaka.txt']

rule: input: TARGETS

path = "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm"
s_studyid = 'GECOPD'

rule reaka_notebook:
    input:
        pop_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.fam',
        chrom_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.bim',
        kin_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.kin',
        sexcheck_file= f'{path}/tmp/{s_studyid}_annotated_plink_merged.sexcheck',
        hwe_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.hwe',
        frq_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.frq',
        het_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.het',
        imiss_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.imiss',
        lmiss_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.lmiss',
        kin0_file = f'{path}/tmp/{s_studyid}_king.kin0'
    output: "tmp/flags_reaka.txt"
    params:
        i_miss_threshold=0.1,
        l_miss_threshold=0.01,
    notebook: "../notebooks/GECOPDNotebook1.ipynb"

rule remra_notebook:
    input:
        pop_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.fam',
        chrom_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.bim',
        kin_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.kin',
        sexcheck_file= f'{path}/tmp/{s_studyid}_annotated_plink_merged.sexcheck',
        hwe_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.hwe',
        frq_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.frq',
        het_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.het',
        imiss_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.imiss',
        lmiss_file = f'{path}/tmp/{s_studyid}_annotated_plink_merged.lmiss',
        kin0_file = f'{path}/tmp/{s_studyid}_king.kin0',
    output: "tmp/flags_remra.txt"
    params:
        i_miss_threshold=0.1,
        l_miss_threshold=0.01,
    notebook: "../notebooks/GECOPDNotebook.ipynb"
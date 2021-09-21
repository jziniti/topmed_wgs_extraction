path = "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm/tmp"
s_studyid = 'GECOPD'

rule reaka_notebook:
    input:
        pop_file = "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm/tmp/sav/GECOPD_annotated_plink_merged.fam",
        chrom_file = "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm/tmp/sav/GECOPD_annotated_plink_merged.bim.original",
        kin_file = f'{path}/{s_studyid}_annotated_plink_merged.kin',
        sexcheck_file= f'{path}/{s_studyid}_annotated_plink_merged.sexcheck',
        hwe_file = f'{path}/{s_studyid}_annotated_plink_merged.hwe',
	frq_file = f'{path}/{s_studyid}_annotated_plink_merged.frq',
	het_file = f'{path}/{s_studyid}_annotated_plink_merged.het',
	imiss_file = f'{path}/{s_studyid}_annotated_plink_merged.imiss',
	lmiss_file = f'{path}/{s_studyid}_annotated_plink_merged.lmiss',
	kin0_file = "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10.cdnm/tmp/GECOPD_king.kin0"
    output: "tmp/reaka_output.txt"
    log:
        notebook="notebook_logs/reaka.ipynb"
    params:
        Error_one=1,
        Error_two=0.5,
        Kinship=0.354,
        SNPSEX=0,
        P=0.000001
	F=0.2
    conda: "../envs/cdnm-jupyter-python-3.7.6.yaml"
    notebook: "../notebooks/GECOPDNotebook1.ipynb"

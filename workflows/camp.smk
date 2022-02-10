module camp_methylation_concordance:
    snakefile: "../850K_qc_pipeline/workflows/external-reference-concordance.smk"
    config: config['CAMP']
    # config: config

use rule * from camp_methylation_concordance as camp_mrc_* 

#use rule write_extract_file from camp_methylation_concordance as camp_mrc_write_extract_file with:
#    output: temp(TMP/'{s_studyid}_extract.txt')

use rule compare_genotypes from camp_methylation_concordance as camp_mfc_compare_genotypes with:
    input:
        control_snp_calls=lambda w: expand(Path(config['CAMP']['methylation_freeze_path'])/f"control_snp_genotypes.{{idx}}.txt", idx=range(config[w.s_studyid]['methylation_batch_count'])),
        ped='/proj/regeps/regep00/studies/CAMP/data/dna/whole_genome/CAMP_FULL_MERGE_UMich/data/freezes/20210928/CAMP_civic_850K_controls_snps.ped',
        map='/proj/regeps/regep00/studies/CAMP/data/dna/whole_genome/CAMP_FULL_MERGE_UMich/data/freezes/20210928/CAMP_civic_850K_controls_snps.map',
        assumed=srcdir("../850K_qc_pipeline/assume.txt"),
        snpmap=srcdir("../850K_qc_pipeline/epic_control_snp_list_GRCh38.p7.txt"),
    output:
        kin=temp(TMP/'{s_studyid}_methylation_gwas_concordance.kin'),
        kin0=temp(TMP/'{s_studyid}_methylation_gwas_concordance.kin0'),

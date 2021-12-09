import pandas as pd

if __name__ == '__main__':
    ### Split the data by race
    df_pheno = pd.read_csv(snakemake.input.pheno, delimiter='\t')
    df_fam = pd.read_csv(snakemake.input.fam, delim_whitespace=True, names=('FID', 'IID', 'DID', 'MID', 'SEX', 'AFF'))
    df_stashq = pd.read_csv(snakemake.input.stashq_file, delim_whitespace=True, names=('ALIASID', '_index_', 'S_SAMPLEID', 'S_SUBJECTID', 'S_STUDYID'))
    df_idmap = pd.read_csv(snakemake.input.idmap_file)


    df_fam_annotated = df_fam.merge(df_stashq, how='left', left_on='IID', right_on='ALIASID')
    df_fam_annotated = df_fam_annotated.merge(df_idmap, how='left', on='S_SUBJECTID')
    df_fam_annotated = df_fam_annotated.merge(df_pheno, how='left', left_on='Pheno Id', right_on='sid')

    if int(snakemake.wildcards.race_code) == 2:
        df_fam_with_matched_race = df_fam_annotated[df_fam_annotated.race == int(snakemake.wildcards.race_code)]
    else:
        df_fam_with_matched_race = df_fam_annotated[df_fam_annotated.race != 2]
    df_fam_with_matched_race.to_csv(snakemake.output.keep, index=False, sep='\t', header=False, columns=('FID', 'IID'))
        
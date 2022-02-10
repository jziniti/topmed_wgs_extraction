import pandas as pd

if __name__ == '__main__':
    fam = pd.read_csv(snakemake.input.fam, delim_whitespace=True, names=('FID', 'IID', 'DID', 'MID', 'SEX', 'AFF'))
    
    df_dict = {'S_SUBJECTID': [],}

    s_studyid = snakemake.wildcards.s_studyid
    if s_studyid == 'ECLPSE':
        s_studyid = 'ECLIPSE'

    with open(snakemake.input.stashq) as fh:
        for line in fh:
            if '=' not in line: continue
            (key, value) = line.strip().split('=')
            if key == 'Subject.S_SUBJECTID':
                s_subjectid = value
            elif key == 'Cohort0.STUDYID' and value == s_studyid:
                df_dict['S_SUBJECTID'].append(s_subjectid)

    study_subset = pd.DataFrame(df_dict)
    study_subset = study_subset.merge(fam, left_on='S_SUBJECTID', right_on='IID', how='left')
    subject_prefix = snakemake.params.get('subject_prefix', '')
    study_subset['prefixed_subject_id'] = study_subset['IID'].apply(lambda w: f'{subject_prefix}{w}')
    study_subset.to_csv(snakemake.output.fam, sep='\t', index=False, header=False, columns=('FID', 'prefixed_subject_id', 'DID', 'MID', 'SEX', 'AFF'))
    # study_subset.to_csv(snakemake.output.fam, sep='\t', index=False, header=False, columns=('FID', 'IID', 'DID', 'MID', 'SEX', 'AFF'))
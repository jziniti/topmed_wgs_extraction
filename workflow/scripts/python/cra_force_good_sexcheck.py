import pandas as pd

if __name__ == '__main__':
    fam = pd.read_csv(snakemake.input.fam, delim_whitespace=True, names=('U_PEDIGREEID', 'S_SUBJECTID', 'DID', 'MID', 'SEX', 'AFF'))
    
    gender_map = {
        'S_SUBJECTID': [],
        'GENDER': [],
    }
    with open(snakemake.input.stashq) as fh:
        for line in fh:
            if '=' not in line: continue
            (key, value) = line.strip().split('=', 1)
            if key == 'Subject.S_SUBJECTID':
                stid = value
            elif key == 'Subject.GENDERFLAG':
                gender = value
                gender_map['S_SUBJECTID'].append(stid)
                gender_map['GENDER'].append(gender)
    gender_df = pd.DataFrame.from_dict(gender_map, orient='columns')
    ## fam['S_SUBJECTID'] = fam['IID'].str.split('_')[1]
    ## fam[['U_PEDIGREEID','S_SUBJECTID']] = fam.IID.str.split('_', expand=True)

    fam = fam.merge(gender_df, on='S_SUBJECTID')
    fam = fam.rename(columns={'S_SUBJECTID': 'SampleLabel', 'GENDER': 'ObservedGender'})
    fam.to_csv(snakemake.output[0], index=False, columns=('SampleLabel', 'ObservedGender'))
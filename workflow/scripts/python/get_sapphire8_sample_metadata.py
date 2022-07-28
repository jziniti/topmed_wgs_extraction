import pandas as pd
from pathlib import Path

from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy.orm import sessionmaker

sapphire8_md = MetaData(schema='SAPPHIRE8')
engine = create_engine('oracle://sapphire_r:sapphire_r0@172.27.104.202/labvprd')
Session = sessionmaker(bind=engine)
session = Session()

SampleFamily = Table('S_SAMPLEFAMILY', sapphire8_md, autoload=True, autoload_with=engine)
Sample = Table('S_SAMPLE', sapphire8_md, autoload=True, autoload_with=engine)

if __name__ == '__main__':
    plate_dfs = []
    for manifest in snakemake.input.plate_manifests:
        plate_name = Path(manifest).stem.split('_')[0]
        plate_df = pd.read_csv(manifest)
        plate_df['PLATEID'] = plate_name
        plate_dfs.append(plate_df)
    plate_metadata = pd.concat(plate_dfs)
    
    all_sample_ids = plate_metadata['SAMPLEID'].unique()

    print(Sample.columns)
    sample_df = pd.read_sql_table(
        's_sample',
        schema='SAPPHIRE8',
        con=engine,
        parse_dates=[],
        columns=['s_sampleid','sampledesc', 'sampletypeid', 'batchid', 'samplefamilyid', 'activeflag']
    )
    sample_df = sample_df[sample_df['s_sampleid'].isin(all_sample_ids)]

    print(SampleFamily.columns)
    sample_family_df = pd.read_sql_table(
        's_samplefamily',
        schema='SAPPHIRE8',
        con=engine,
    )
    
    annotated_sample_df = sample_df.merge(sample_family_df, left_on='samplefamilyid', right_on='s_samplefamilyid', how='left', suffixes=(None, '_samplefamily'))
    #print(annotated_sample_df)
    annotated_sample_df.to_csv(snakemake.output[0], index=False)

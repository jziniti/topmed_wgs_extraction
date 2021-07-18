import pandas as pd


#### BOILERPLATE
from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy.orm import sessionmaker
sapphire8_md = MetaData(schema='SAPPHIRE8')
engine = create_engine('oracle://sapphire_r:sapphire_r0@172.27.104.202:1521/labvprd')
Session = sessionmaker(bind=engine)
session = Session()
SDIAlias = Table('SDIALIAS', sapphire8_md, autoload=True, autoload_with=engine)
##### BOILERPLATE

alias_prefix = snakemake.params.get('alias_prefix', 'NWD')
s_studyid = snakemake.params.get('s_studyid', None)
if s_studyid == 'ALL': s_studyid = None
if s_studyid:
    nwdids = pd.read_sql(session.query(SDIAlias).where(SDIAlias.c.aliasid.like(f'{alias_prefix}%')).filter(SDIAlias.c.keyid2==s_studyid).statement, engine)
else:
    nwdids = pd.read_sql(session.query(SDIAlias).where(SDIAlias.c.aliasid.like(f'{alias_prefix}%')).statement, engine)
## nwdids = nwdids.rename(columns={'keyid1':'S_SAMPLEID'})
nwdids.to_csv(snakemake.output.nwdids, index=False, columns=['aliasid',])

import cx_Oracle
from AIMS.SAPPHIRE import SDIAlias, Sample, SampleFamily
from AIMS.Generations import util

if __name__ == '__main__':
    cnx = util.connect('sapphire_r@labvprd')

    ## nwd_aliases = SDIAlias.
    aliases = util.doSelect("SELECT * FROM SAPPHIRE8.SDIALIAS WHERE ALIASID LIKE 'TOR%'", cnx)
    print('ALIAS,SAMPLEID,STUDYID,SUBJECTID')
    for alias in aliases:
        aid = alias['ALIASID']
        sid = alias['KEYID1']
        study = alias['KEYID2']


        if not sid.startswith('S-'): continue ## Some samples are mistakenly mapped to ST-ID.  ignore these ...
        sample = Sample()
        sample.select(cnx, S_SAMPLEID=sid)

        sf = SampleFamily()
        sf.select(cnx, S_SAMPLEFAMILYID=sample.SAMPLEFAMILYID)

        stid = sf['SUBJECTID']

        print('%s,%s,%s,%s' % (aid, sid, study, stid))
        

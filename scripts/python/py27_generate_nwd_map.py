import cx_Oracle
from AIMS.SAPPHIRE import SDIAlias, Sample, SampleFamily
from AIMS.Generations import util

if __name__ == '__main__':
    cnx = util.connect('sapphire_r@labvprd')

    ## nwd_aliases = SDIAlias.
    nwd_aliases = util.doSelect("SELECT * FROM SAPPHIRE8.SDIALIAS WHERE ALIASID LIKE 'NWD%'", cnx)
    print('NWDID,SAMPLEID,STUDYID,SUBJECTID')
    for alias in nwd_aliases:
        nwdid = alias['ALIASID']
        sid = alias['KEYID1']
        study = alias['KEYID2']

        sample = Sample()
        sample.select(cnx, S_SAMPLEID=sid)

        sf = SampleFamily()
        sf.select(cnx, S_SAMPLEFAMILYID=sample.SAMPLEFAMILYID)

        stid = sf['SUBJECTID']

        print('%s,%s,%s,%s' % (nwdid, sid, study, stid))
        

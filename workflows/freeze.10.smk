STUDIES = ['CAMP',
           'CRA',
           'ECLPSE',
           'EOCOPD',
           'GECOPD',
           'GLAXO',
           'LTCOPD',
           'LTRC',
           'PLCOPD'
]
TARGETS = expand("tmp/{s_studyid}.nwds.txt", s_studyid=STUDIES)
TARGETS.append('tmp/ALL.nwds.txt')

rule: input: TARGETS

rule sapphire_nwd_list:
    input: 
    output: nwdids=temp("tmp/{s_studyid}.nwds.txt")
    conda: "../envs/sapphire8.yaml"
    params:
        alias_prefix='NWD',
        s_studyid='{s_studyid}',
    script: "../scripts/python/read_sapphire_aliases.py"


from snakemake.utils import min_version
min_version("6.0")

TMP = Path('tmp')

TARGETS.append('ANNOTATED_MANIFEST.csv')
rule: input: TARGETS

module multiomics: 
    snakefile: "../topmed_multiomics_qc/workflows/multiomics_qc.smk" ### FIXME: snakefile: "https://URL/ would be good here?"
    config: config['GECOPD']
    
use rule * from multiomics

use rule add_manual_flags from multiomics with:
    output:
        #df='multiomics/{s_studyid}/ANNOTATED_MANIFEST.csv',
        df='{tmpdir}/ANNOTATED_MANIFEST.csv',

use rule read_gender_observations from multiomics with:
    input:
        lambda w: config[w.s_studyid]['datasets'][int(w.dataset)]['gender_observations'],
    output:
        df=temp("multiomics/{s_studyid}/{dataset}.gender_observations"),

use rule read_concordance_data from multiomics with:
    input:
        lambda w: config[w.s_studyid]['concordance_observations'][int(w.index)]['path'],
    output:
        df=temp("multiomics/{s_studyid}/{index}.concordance_data"),

use rule read_dataset_manifest from multiomics with:
    input:
        lambda w: config[w.s_studyid]['datasets'][int(w.dataset)]['manifest_path'],
    output:
        cdnm_sids=temp("multiomics/{s_studyid}/{dataset}.cdnm_sids"),
        cdnm_aliases=temp("multiomics/{s_studyid}/{dataset}.cdnm_aliases"),
        manifest=temp("multiomics/{s_studyid}/{dataset}.manifest"), 

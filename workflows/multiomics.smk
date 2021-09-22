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
        df='multiomics/{s_studyid}/ANNOTATED_MANIFEST.csv',
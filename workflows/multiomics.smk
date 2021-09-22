from snakemake.utils import min_version
min_version("6.0")

TMP = Path('tmp')

TARGETS = ['ANNOTATED_MANIFEST.csv',]
rule: input: TARGETS

module multiomics:
    snakefile: "../topmed_multiomics_qc/workflows/multiomics_qc.smk"
    config: config['GECOPD']

use rule * from multiomics
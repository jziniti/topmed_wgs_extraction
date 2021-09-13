from snakemake.utils import min_version
min_version("6.0")

TMP = Path('tmp')

rule: input: [TMP/'fake.kin0']

module methylation_ref_concordance:
    snakefile: "../850K_qc_pipeline/workflows/external-reference-concordance.smk"
    config: config

use rule * from methylation_ref_concordance as mrc_* 

use rule compare_genotypes from methylation_ref_concordance as mrc_compare_genotypes with:
    input:
        control_snp_calls=expand("/proj/regeps/regep00/studies/GECOPD/data/methylation/TopMed/data/freezes/20200706/LEVEL2/control_snp_genotypes.{idx}.txt", idx=range(config['batch_count'])),


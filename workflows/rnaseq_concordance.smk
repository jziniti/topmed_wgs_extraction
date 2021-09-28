from snakemake.utils import min_version
min_version("6.0")

try: TMP
except NameError: TMP = Path('tmp/')

try: TARGETS
except NameError: TARGETS = []

try: STUDIES
except NameError: STUDIES = sorted(config.keys())

KING_DIR = TMP/'king'

#for s_studyid in STUDIES:
#    if config[s_studyid].get('run_rna_concordance', False):
#        TARGETS.append(TMP/f'{s_studyid}_rna_king_results_summary.html')
#        TARGETS.append(TMP/f'{s_studyid}_rna_king_results_summary.csv')

rule: input: TARGETS

module rnaseq_concordance:
    snakefile: "../rnaseq_snp_concordance/workflows/concordance_only.smk"
    config: config['GECOPD']

use rule * from rnaseq_concordance as rna_* 

use rule generate_king_report from rnaseq_concordance as rna_generate_king_report with:
    output:
        html_file=TMP/"{s_studyid}_rna_king_results_summary.html",
        csv_file=TMP/"{s_studyid}_rna_king_results_summary.csv",

use rule bcftools_setid from rnaseq_concordance with:
     input: vcf=lambda w: f'/proj/regeps/regep00/studies/COPDGene/analyses/rejpz/bcbio_ea/data/bcbio/{w.bambase}/final/{w.bambase}/{w.bambase}-freebayes.vcf.gz',
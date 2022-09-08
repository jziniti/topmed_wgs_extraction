from snakemake.utils import min_version
min_version("6.0")

TMP = Path('tmp')

# TARGETS.append('ANNOTATED_MANIFEST.csv')

config['input_files'] = [TMP/"{s_studyid}_annotated_plink_merged.vcf.gz",]

module cdnm_pca_pipeline:
    snakefile: "../cdnm_pca_pipeline/workflows/pre-process.smk"
    config: config

use rule * from cdnm_pca_pipeline

use rule filter_input_genotyped_maf from cdnm_pca_pipeline with:
    input:
        vcf=TMP/"{s_studyid}_annotated_plink_merged.vcf.gz"
    output:
        vcf=temp(TMP/"{s_studyid}_filtered.vcf.gz")

use rule concat_passthrough from cdnm_pca_pipeline with:
        input:
            vcf=TMP/"{s_studyid}_filtered.vcf.gz"
        output:
            vcf=TMP/"{s_studyid}_laser_trace_input.vcf.gz"

use rule manual_run_on_umich from cdnm_pca_pipeline with:
    input:
        laser_input=TMP/"{s_studyid}_laser_trace_input.vcf.gz"
    output:
        touch(TMP/"{s_studyid}_ready_for_umich.done")


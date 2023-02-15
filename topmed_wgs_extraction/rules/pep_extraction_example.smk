import pathlib

pepfile: config['pepfile']
RAW_WGS_BASE_PATH = pathlib.Path(pep.config['base_path'])
RAW_WGS_FILES = pep.config['files']
SAMPLE_ID_STRING = ','.join([sample['sample_name'] for sample in pep.samples if sample['action'] == 'keep'])

TMP = pathlib.Path(config.get('scratch_dir', 'tmp'))
OUT = pathlib.Path(config.get('extract_dir', '.')).absolute()
LOG_DIR = pathlib.Path(config.get('logdir', './logs/'))

PAR1_RANGE = 'chrX:10001-2781479'
PAR2_RANGE = 'chrX:155701383-156030895'

### A default set of targets for Snakemake
ALL = []
for filename in RAW_WGS_FILES:
    ALL.append(str(OUT/filename))
    if 'chrX' in filename:
        ALL.append(OUT/filename.replace('chrX', 'chrXY'))

rule all:
    input: ALL

rule extract:
    input: bcf=RAW_WGS_BASE_PATH/"{bcf_name}"
    output: bcf=OUT/"{bcf_name}"
    conda: "../envs/bcftools.yaml"
    log: LOG_DIR/'extract.{bcf_name}.log'
    params: samples=SAMPLE_ID_STRING
    shell: "bcftools view -s {params.samples} -i 'FILTER=\"PASS\"' -c 1 --force-samples {input.bcf} \
            | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O b -o {output.bcf}"

rule extract_X:
    input: bcf=RAW_WGS_BASE_PATH/"freeze.10a.chrX.pass_and_fail.gtonly.minDP10.bcf"
    output: bcf=OUT/"freeze.10a.chrX.pass_and_fail.gtonly.minDP10.bcf"
    conda: "../envs/bcftools.yaml"
    log: LOG_DIR/'extract.freeze.10a.chrX.pass_and_fail.gtonly.minDP10.bcf.log'
    params:
        samples=SAMPLE_ID_STRING,
        par1_range=PAR1_RANGE,
        par2_range=PAR2_RANGE,
    shell: """bcftools view {input.bcf} \
                          --targets ^{params.par1_range},{params.par2_range} \
                          -s {params.samples} \
                          -i 'FILTER=\"PASS\"' \
                          -c 1 \
                          --force-samples \
              | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O b -o {output.bcf}"""

rule extract_XY:
    input: bcf=RAW_WGS_BASE_PATH/"freeze.10a.chrX.pass_and_fail.gtonly.minDP10.bcf"
    output: bcf=OUT/"freeze.10a.chrXY.pass_and_fail.gtonly.minDP10.bcf"
    conda: "../envs/bcftools.yaml"
    log: LOG_DIR/'extract.freeze.10a.chrXY.pass_and_fail.gtonly.minDP10.bcf.log'
    params:
        samples=SAMPLE_ID_STRING,
        par1_range=PAR1_RANGE,
        par2_range=PAR2_RANGE,
    shell: """bcftools view {input.bcf} \
                          --regions {params.par1_range},{params.par2_range} \
                          -s {params.samples} \
                          -i 'FILTER=\"PASS\"' \
                          --force-samples \
                          -c 1 \
                          -O b \
                          -o {output.bcf}"""


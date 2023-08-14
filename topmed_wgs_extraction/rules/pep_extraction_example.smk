import pathlib

pepfile: config['pepfile']
SAMPLE_TABLE_PATH = config['sample_table']
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
    ALL.append(str(OUT/f'{filename}.csi'))
    if 'chrX' in filename:
        ALL.append(OUT/filename.replace('chrX', 'chrXY'))
        ALL.append(OUT/f"{filename.replace('chrX', 'chrXY')}.csi")

wildcard_constraints:
    bcf_name='\S+.bcf'
    
rule all:
    input: ALL

rule index:
    input: bcf=OUT/"{bcf_name}"
    output: csi=OUT/"{bcf_name}.csi"
    conda: "../envs/bcftools.yaml"
    shell: "bcftools index {input.bcf} --force --csi"

rule extract:
    input: bcf=RAW_WGS_BASE_PATH/"{bcf_name}"
    output: bcf=OUT/"{bcf_name,\w+.bcf}"
    conda: "../envs/bcftools.yaml"
    log: LOG_DIR/'extract.{bcf_name}.log'
    params: samples=SAMPLE_ID_STRING
    shell: "bcftools view -s {params.samples} -i 'FILTER=\"PASS\"' -c 1 --force-samples {input.bcf} \
            | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -O b -o {output.bcf}"

rule extract_X:
    input:
         bcf=RAW_WGS_BASE_PATH/"freeze.10a.chrX.pass_and_fail.gtonly.minDP10.bcf",
         sample_table=SAMPLE_TABLE_PATH,
    output: bcf=OUT/"freeze.10a.chrX.pass_and_fail.gtonly.minDP10.bcf"
    conda: "../envs/bcftools+python.yaml"
    log: LOG_DIR/'extract.freeze.10a.chrX.pass_and_fail.gtonly.minDP10.bcf.log'
    params:
        samples=SAMPLE_ID_STRING,
        par1_range=PAR1_RANGE,
        par2_range=PAR2_RANGE,
        het_hap_filter_script=srcdir('../scripts/python/chrX_remove_het_male_snps.py'),
        # sample_table=pep.sample_table,
    shell: """bcftools view {input.bcf} \
                          --targets ^{params.par1_range},{params.par2_range} \
                          -s {params.samples} \
                          -i 'FILTER=\"PASS\"' \
                          -c 1 \
                          --force-samples \
              | python {params.het_hap_filter_script} {input.sample_table} | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -O b -o {output.bcf}"""

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
              | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -O b -o {output.bcf}"""


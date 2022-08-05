import pathlib

### Everything we know about the dataset should come from the PEP. If
### not in there, and you think it should be, ask someone upstream of
### you to add it
pepfile: config['pepfile']
RAW_WGS_BASE_PATH = pathlib.Path(pep.config['base_path'])
RAW_WGS_FILES = pep.config['files']
SAMPLE_ID_STRING = ','.join([sample['sample_name'] for sample in pep.samples if sample['action'] == 'keep'])

OUT = pathlib.Path(config.get('extract_dir', '.')).absolute()

### A default set of targets for Snakemake
ALL = []
for filename in RAW_WGS_FILES:
    ALL.append(str(OUT/filename))

rule all:
    input: ALL
                      
rule extract:
    input:
        bcf=RAW_WGS_BASE_PATH/"{bcf_name}"
    output:
        bcf=OUT/"{bcf_name}"
    conda: "../envs/bcftools.yaml"
    params: samples=SAMPLE_ID_STRING
    # shell: "bcftools view -s {params.samples} -i 'FILTER=\"PASS\"' -c 1 -O b -m2 -M2 --force-samples --types snps {input.bcf} -o {output.bcf}"
    shell: "bcftools view -s {params.samples} -i 'FILTER=\"PASS\"' -O b --force-samples {input.bcf} -o {output.bcf}"

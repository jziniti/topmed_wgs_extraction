import pathlib

LOWQUAL_BASE_DIR = pathlib.Path('/d/tmp2/regeps/regep00/studies/COPDGene/data/freezes/lowqual/')
LOWQUAL_BATCHES = ['2018.0819.weiss.rejects.59', '2018.0819.weiss.rejects.60']

TMP = Path('tmp')

TARGETS = [TMP/'NWD999786']

rule: input: TARGETS

rule extract:
    input: cram=LOWQUAL_BASE_DIR/LOWQUAL_BATCHES[0]/'{ALIAS}.src.cram'
    output: TMP/'{ALIAS}'
    conda: "../envs/samtools.yaml"
    shell: "samtools query -l {input.cram}"
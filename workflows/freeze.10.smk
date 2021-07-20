STUDIES = config['studies']
TARGETS = expand("{s_studyid}.bed", s_studyid=STUDIES)

rule done: input: TARGETS

rule read_sample_annotation:
    input: 
    output: nwdids=temp("tmp/{s_studyid}.dcc.txt")
    conda: "../envs/sapphire8.yaml"
    script: "../scripts/python/read_sapphire_aliases.py"

rule sapphire_nwd_list:
    input: 
    output: nwdids=temp("tmp/{s_studyid}.sapphire.txt")
    conda: "../envs/sapphire8.yaml"
    params:
        alias_prefix='NWD',
        s_studyid='{s_studyid}',
    script: "../scripts/python/read_sapphire_aliases.py"

rule sample_sanity_check:
    input:
        "tmp/{s_studyid}.dcc.txt",
        "tmp/{s_studyid}.sapphire.txt",
    output: nwdids=temp("tmp/{s_studyid}.nwds.txt")
    conda: "../envs/sapphire8.yaml"
    shell: "echo"

rule extract:
    input: "tmp/{s_studyid}.nwds.txt"
    output: temp("tmp/{s_studyid}.extract.vcf.gz")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools"

rule stashq:
    input: "tmp/{s_studyid}.nwds.txt"
    output: temp("tmp/{s_studyid}.stashq")
    conda: "../envs/stashq.yaml"
    shell: "stashq alias `cat {input}` > {output}"

rule annotate:
    input:
        stashq="tmp/{s_studyid}.stashq",
        vcf="tmp/{s_studyid}.extract.vcf.gz",
    output: temp("tmp/{s_studyid}.annotate.vcf.gz")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools"

rule variant_qc:
    input: "tmp/{s_studyid}.annotate.vcf.gz"
    output: temp("tmp/{s_studyid}.variants.keep")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools"

rule sample_qc:
    input: "tmp/{s_studyid}.annotate.vcf.gz"
    output: temp("tmp/{s_studyid}.samples.keep")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools"

rule filter:
    input:
        "tmp/{s_studyid}.annotate.vcf.gz",
        "tmp/{s_studyid}.samples.keep",
        "tmp/{s_studyid}.variants.keep",
    output: temp("{s_studyid}.bed")
    conda: "../envs/bcftools.yaml"
    shell: "bcftools"

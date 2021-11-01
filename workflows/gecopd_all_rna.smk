TARGETS = ['bam_sample_map.tsv']

BAMS_BY_BASE = {}
INPUT_BAMS = {}

rule: input: TARGETS 

rule list_bam_paths:
    input: "bam_paths.txt"
    output: temp("tmp/all_bam_files.txt")
    shell: "ls -1 `cat {input}` > {output}"

rule md5sums:
    input: "tmp/all_bam_files.txt"
    output: temp("tmp/all_bam_files.md5sums")
    shell: "md5sum `cat {input}` > {output}"

rule prepare_bam_sample_map:
    input: "tmp/all_bam_files.md5sums"
    output: "bam_sample_map.tsv"
    script: "../scripts/python/prepare_bam_sample_map.py"


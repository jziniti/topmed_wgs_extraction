
rule cra_ref_stashq_subject:
    input: fam="/proj/regeps/regep00/studies/CRA/data/dna/whole_genome/2013_CRA_GWAS_UMich/data/freezes/20190523_CIVIC/2013_CRA_civic.fam",
    output: temp(TMP/"2013_CRA_civic.stashq.txt")
    conda: "../cdnm/envs/stashq.yaml"
    shell: "stashq subject `awk '{{print $2}}' {input.fam} | cut -d \"_\" -f 2` > {output}"

rule force_good_sexcheck:
    input:
        fam="/proj/regeps/regep00/studies/CRA/data/dna/whole_genome/2013_CRA_GWAS_UMich/data/freezes/20190523_CIVIC/2013_CRA_civic.fam",
        stashq="tmp/2013_CRA_civic.stashq.txt"
    output: "tmp/2013_CRA_civic.forced_good_sexcheck"
    script: "../scripts/python/cra_force_good_sexcheck.py" 
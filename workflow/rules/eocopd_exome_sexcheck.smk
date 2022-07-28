
rule eocopd_convert_bim_chrpos:
    input:
        bim="/proj/regeps/regep00/studies/EOCOPD/data/dna/exome_chip/BWH_Silverman_Exome6k/data/freezes/20150315/exome6kSubAndMarkCleanV03.bim",
    output:
        bim="tmp/exome6kSubAndMarkCleanV03_chrpos.bim"
    shell: "awk '{{print $1,$1\":\"$4,$3,$4,$5,$6}}' {input.bim} > {output.bim}"

rule eocopd_exome_prep:
    input:
        bed='/proj/regeps/regep00/studies/EOCOPD/data/dna/exome_chip/BWH_Silverman_Exome6k/data/freezes/20150315/exome6kSubAndMarkCleanV03.bed',
        bim='tmp/exome6kSubAndMarkCleanV03_chrpos.bim',
        fam='/proj/regeps/regep00/studies/EOCOPD/data/dna/exome_chip/BWH_Silverman_Exome6k/data/freezes/20150315/exome6kSubAndMarkCleanV03.fam',
    output:
        bed='tmp/exome6kSubAndMarkCleanV03.bed',
        bim='tmp/exome6kSubAndMarkCleanV03.bim',
        fam='tmp/exome6kSubAndMarkCleanV03.fam',
    conda: "../cdnm/envs/plink.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out tmp/exome6kSubAndMarkCleanV03 --make-bed"
        
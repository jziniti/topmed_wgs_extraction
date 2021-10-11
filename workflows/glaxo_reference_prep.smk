rule glaxo_convert_bim_chrpos:
    input:
        bim="/proj/regeps/regep00/studies/ICGN/data/dna/whole_genome/AxiomGenotyping/data/freezes/20191011/ICGN_AxiomGenotyping.bim",
    output:
        bim="tmp/ICGN_AxiomGenotyping_munged.bim"
    shell: "awk '{{print $1,$1\":\"$4,$3,$4,$5,$6}}' {input.bim} > {output.bim}"

rule prepare_reference_for_concordance:
    input:
        bed="/proj/regeps/regep00/studies/ICGN/data/dna/whole_genome/AxiomGenotyping/data/freezes/20191011/ICGN_AxiomGenotyping.bed",
        bim="tmp/ICGN_AxiomGenotyping_munged.bim",
        fam="/proj/regeps/regep00/studies/ICGN/data/dna/whole_genome/AxiomGenotyping/data/freezes/20191011/ICGN_AxiomGenotyping.fam",
    output:
        bed="tmp/ICGN_AxiomGenotyping_chrpos.bed",
        bim="tmp/ICGN_AxiomGenotyping_chrpos.bim",
        fam="tmp/ICGN_AxiomGenotyping_chrpos.fam",
    params: out="tmp/ICGN_AxiomGenotyping_chrpos"
    conda: "../cdnm/envs/plink.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out {params.out} --make-bed"
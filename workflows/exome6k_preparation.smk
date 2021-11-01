rule exome6k_study_subset:
    input:
        fam="/proj/regeps/regep00/studies/EOCOPD/data/dna/exome_chip/BWH_Silverman_Exome6k/data/freezes/20150315/exome6kSubAndMarkCleanV03.fam",
        stashq="tmp/exome6k.stashq.txt",
    output:
        fam="tmp/{s_studyid}_exome6kSubAndMarkCleanV03.fam",
    script: "../scripts/python/exome6k_study_subset.py"

rule axiom_study_subset:
    input:
        fam="/proj/regeps/regep00/studies/ICGN/data/dna/whole_genome/AxiomGenotyping/data/freezes/20191011/ICGN_AxiomGenotyping.fam",
        stashq="tmp/Axiom.stashq.txt",
    output:
        fam="tmp/{s_studyid}_ICGN_AxiomGenotyping.fam",
    params:
        subject_prefix="ax:"
    script: "../scripts/python/exome6k_study_subset.py"

rule axiom_modify_fam_ax:
    input: fam="tmp/ICGN_AxiomGenotyping_chrpos.fam"
    output: fam="tmp/ICGN_AxiomGenotyping_w_ax_chrpos.fam"
    shell: "cat {input.fam} | awk '{{print $1,\"ax:\"$2,$3,$4,$5,$6}}' > {output.fam}"

rule axiom_recode_fam_ax:
    input:
        bed="tmp/ICGN_AxiomGenotyping_chrpos.bed",
        bim="tmp/ICGN_AxiomGenotyping_chrpos.bim",
        fam="tmp/ICGN_AxiomGenotyping_w_ax_chrpos.fam",
    output:
        bed="tmp/ICGN_AxiomGenotyping_ax_chrpos.bed",
        bim="tmp/ICGN_AxiomGenotyping_ax_chrpos.bim",
        fam="tmp/ICGN_AxiomGenotyping_ax_chrpos.fam",
    params: out="tmp/ICGN_AxiomGenotyping_ax_chrpos"
    conda: "../cdnm/envs/plink.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out {params.out} --make-bed"

use rule duplicate_vs_reference from king as axiom_reference_concordance with:
    input:
        rbed="tmp/ICGN_AxiomGenotyping_ax_chrpos.bed",
        rbim="tmp/ICGN_AxiomGenotyping_ax_chrpos.bim",
        rfam="tmp/ICGN_AxiomGenotyping_ax_chrpos.fam",
        qbed="tmp/{s_studyid}_annotated_plink_merged_refformat.bed",
        qbim="tmp/{s_studyid}_annotated_plink_merged_refformat.bim",
        qfam="tmp/{s_studyid}_annotated_plink_merged_refformat.fam",
    output:
        kin=TMP/"{s_studyid}_axiom_concordance.con",
    resources:
        mem_free="200G",
        mem_mb="2000",
    params:
        Q=lambda w: f"tmp/{w.s_studyid}_annotated_plink_merged_refformat",
        R="tmp/ICGN_AxiomGenotyping_ax_chrpos",
        prefix=lambda w: TMP/f"{w.s_studyid}_axiom_concordance",

use rule duplicate_vs_reference from king as exome6k_reference_concordance with:
    input:
        rbed="tmp/exome6kSubAndMarkCleanV03_chrpos.bed",
        rbim="tmp/exome6kSubAndMarkCleanV03_chrpos.bim",
        rfam="tmp/exome6kSubAndMarkCleanV03_chrpos.fam",
        qbed="tmp/{s_studyid}_annotated_plink_merged_refformat.bed",
        qbim="tmp/{s_studyid}_annotated_plink_merged_refformat.bim",
        qfam="tmp/{s_studyid}_annotated_plink_merged_refformat.fam",
    output:
        kin=TMP/"{s_studyid}_exome6k_concordance.con",
    resources:
        mem_free="200G",
        mem_mb="2000",
    params:
        Q=lambda w: f"tmp/{w.s_studyid}_annotated_plink_merged_refformat",
        R="tmp/exome6kSubAndMarkCleanV03_chrpos",
        prefix=lambda w: TMP/f"{w.s_studyid}_exome6k_concordance",

"""
## fam="/proj/regeps/regep00/studies/ICGN/data/dna/whole_genome/AxiomGenotyping/data/freezes/20191011/ICGN_AxiomGenotyping.fam",
### stashq subject `cut -d " " -f 2 /proj/regeps/regep00/studies/EOCOPD/data/dna/exome_chip/BWH_Silverman_Exome6k/data/freezes/20150315/exome6kSubAndMarkCleanV03.fam` > ../tmp/exome6k.stashq.txt
### stashq subject `cut -d " " -f 2 /proj/regeps/regep00/studies/ICGN/data/dna/whole_genome/AxiomGenotyping/data/freezes/20191011/ICGN_AxiomGenotyping.fam` > ../tmp/Axiom.stashq.txt
"""
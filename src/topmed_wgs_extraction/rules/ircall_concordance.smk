REFERENCE_DATASET_LOOKUP = {
    'exome': TMP/"exome6kSubAndMarkCleanV03.bed",
    'axiom': TMP/'ICGN_AxiomGenotyping_chrpos.bed',
    'GECOPD': '/proj/regeps/regep00/studies/COPDGene/data/genotype/data/freezes/king_ref_GRCh38/king_ref_GRCh38.bed',
    'CAMP': '/proj/regeps/regep00/studies/CAMP/data/dna/whole_genome/CAMP_FULL_MERGE_UMich/data/freezes/20210928/CAMP_civic.bed',
    'CRA': '/proj/regeps/regep00/studies/CRA/data/dna/whole_genome/2013_CRA_GWAS_UMich/data/freezes/20190523_CIVIC/2013_CRA_civic.bed',
    'ECLPSE': '/proj/regeps/regep00/studies/ECLIPSE/data/imputed/ECLIPSE_HRC1-1/data/freezes/20210929/eclipseHrc11.bed',
}

REFERENCE_FORMAT_LOOKUP = {
    'CAMP': '_civic',
    'CRA': '_civic',
    'ECLPSE': '_civic',
    'GECOPD': '_chrpos',
    'axiom': '_chrpos',
    'exome': '_chrpos',
}

rule ircall_concordance:
    input:
        qbed=lambda w: TMP/f"IRCALL_annotated_plink_merged{REFERENCE_FORMAT_LOOKUP.get(w.reference_label, '')}.bed",
        rbed=lambda w: REFERENCE_DATASET_LOOKUP[w.reference_label]
    output:
        TMP/'IRCALL_{reference_label}_duplicate.con',
    params:
        prefix=lambda w: TMP/f'IRCALL_{w.reference_label}_duplicate',
    conda: "../cdnm/envs/king.yaml"
    shell: "king -b {input.rbed},{input.qbed} --prefix {params.prefix} --duplicate"
    #wrapper: "https://raw.githubusercontent.com/CDNMBioinformatics/KING_Wrappers/main/duplicate/2bed/base"

rule axiom_subsample:
    input:
        bed=TMP/f"IRCALL_annotated_plink_merged_chrpos.bed",
        bim=TMP/f"IRCALL_annotated_plink_merged_chrpos.bim",
        fam=TMP/f"IRCALL_annotated_plink_merged_chrpos.fam",
    output:
        bed=TMP/f"IRCALL_annotated_plink_merged_subsample_chrpos.bed",
        bim=TMP/f"IRCALL_annotated_plink_merged_subsample_chrpos.bim",
        fam=TMP/f"IRCALL_annotated_plink_merged_subsample_chrpos.fam",
    conda: "../cdnm/envs/plink.yaml"
    shell: "plink --thin 0.1 --bed {input.bed} --bim {input.bim} --fam {input.fam} --out tmp/IRCALL_annotated_plink_merged_subsample_chrpos --make-bed"

rule ircall_axiom_concordance:
    input:
        qbed=TMP/f"IRCALL_annotated_plink_merged_subsample_chrpos.bed",
        rbed=TMP/'ICGN_AxiomGenotyping_chrpos.bed'
    output:
        TMP/'IRCALL_axiom1_duplicate.con',
    params:
        prefix=TMP/f'IRCALL_axiom1_duplicate',
    conda: "../cdnm/envs/king.yaml"
    shell: "king -b {input.rbed},{input.qbed} --prefix {params.prefix} --duplicate"
    #wrapper: "https://raw.githubusercontent.com/CDNMBioinformatics/KING_Wrappers/main/duplicate/2bed/base"

rule ircall_convert_bim_to_civic:
    input: bim=TMP/"IRCALL_annotated_plink_merged.bim",
    output: bim=temp(TMP/"IRCALL_annotated_plink_merged_c1.bim")
    shell: "awk '{{print $1,$1\":\"$4\":\"$5\":\"$6,$3,$4,$5,$6}}' {input.bim} > {output.bim}"


rule ircall_convert_bim_to_chrpos:
    input: bim=TMP/"IRCALL_annotated_plink_merged.bim"
    output: bim=temp(TMP/"IRCALL_annotated_plink_merged_c2.bim")
    shell: "awk '{{print $1,$1\":\"$4,$3,$4,$5,$6}}' {input.bim} > {output.bim}"

rule ircall_recode_chrpos:
    input:
        bed="tmp/IRCALL_annotated_plink_merged.bed",
        bim="tmp/IRCALL_annotated_plink_merged_c2.bim",
        fam="tmp/IRCALL_annotated_plink_merged.fam",
    output:
        bed="tmp/IRCALL_annotated_plink_merged_chrpos.bed",
        bim="tmp/IRCALL_annotated_plink_merged_chrpos.bim",
        fam="tmp/IRCALL_annotated_plink_merged_chrpos.fam",
    params: out="tmp/IRCALL_annotated_plink_merged_chrpos"
    conda: "../cdnm/envs/plink.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out {params.out} --make-bed"

rule ircall_recode_civic:
    input:
        bed="tmp/IRCALL_annotated_plink_merged.bed",
        bim="tmp/IRCALL_annotated_plink_merged_c1.bim",
        fam="tmp/IRCALL_annotated_plink_merged.fam",
    output:
        bed="tmp/IRCALL_annotated_plink_merged_civic.bed",
        bim="tmp/IRCALL_annotated_plink_merged_civic.bim",
        fam="tmp/IRCALL_annotated_plink_merged_civic.fam",
    params: out="tmp/IRCALL_annotated_plink_merged_civic"
    conda: "../cdnm/envs/plink.yaml"
    shell: "plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --out {params.out} --make-bed"
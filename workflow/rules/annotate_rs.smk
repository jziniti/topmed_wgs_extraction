include: "/udd/rejpz/projects/cdnm/rules/tabix.smk"

TARGETS = ["tmp/GECOPD_annotated_plink_merged_rs_GRCh38.vcf.gz", "tmp/GECOPD_annotated_plink_merged_rs_GRCh38.vcf.gz.tbi"]

rule: input: TARGETS

rule annotate_rs:
    input:
        vcf="tmp/GECOPD_annotated_plink_merged.vcf.gz",
        tbi="tmp/GECOPD_annotated_plink_merged.vcf.gz.tbi",
        # dbsnp="/udd/repsa/regeps/dbsnp/All_20170403.chr1.norm.civicized2.vcf.gz",
        dbsnp="/proj/rerefs/reref00/dbsnp/snp151/GRCh38.p7/VCF/All_20180418.vcf.gz",
    output: vcf="tmp/GECOPD_annotated_plink_merged_rs.vcf.gz"
    conda: "../envs/bcftools.yaml"
    shell: "bcftools annotate {input.vcf} -c CHROM,FROM,TO,ID -a {input.dbsnp} -Oz -o {output.vcf}"

rule chromsome_map:
    output: "tmp/GRCH38_chromosome_map.tsv"
    run:
        with open(output[0], 'w') as fh:
            for chrom in range(1, 26):
                fh.write(f'{chrom}\tchr{chrom}\n')

rule rename_chrs:
    input:
        vcf="tmp/GECOPD_annotated_plink_merged_rs.vcf.gz",
        tbi="tmp/GECOPD_annotated_plink_merged_rs.vcf.gz.tbi",
        chromosome_map="tmp/GRCH38_chromosome_map.tsv",
    output:
        vcf="tmp/GECOPD_annotated_plink_merged_rs_GRCh38.vcf.gz",
    conda: "../envs/bcftools.yaml"
    shell: "bcftools annotate --rename-chrs {input.chromosome_map} {input.vcf} -Oz -o {output.vcf}"

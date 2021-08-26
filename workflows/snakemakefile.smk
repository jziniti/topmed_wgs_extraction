CHRS = list(("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"))



rule done:
	input:
		"king_duplicate.con",
		"sexcheck.sexcheck"

rule extract_bcf:
	input:
		samples="samples",
		bcf="/proj/edith/regeps/regep00/studies/COPDGene/data/wgs/TopMed/data/freezes/freeze.10a/minDP10/freeze.10a.chr{chr}.pass_and_fail.gtonly.minDP10.bcf"
	output:
		out_bcf=temp("CRA_chr{chr}_PASS.bcf")
	shell:
		"""
		bcftools view -S {input.samples} -i 'FILTER=\"PASS\"' -c 1 -O b -m2 -M2 --force-samples --types snps {input.bcf} > {output.out_bcf}
		"""


rule annotate_bcf:
	input:
		"CRA_chr{chr}_PASS.bcf"
	output:
		temp("CRA_annotated_chr{chr}.bcf")
	shell:
		"""
		bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O b -o {output} {input}
		"""

rule transform_plink:
	input:
		"CRA_annotated_chr{chr}.bcf"
	output:
		bed="CRA_annotated_plink_chr{chr}.bed",
		bim="CRA_annotated_plink_chr{chr}.bim",
		fam="CRA_annotated_plink_chr{chr}.fam"
	params:
		chrom='{chr}'
	shell:
		"""
		plink2_v2.00a3LM --bcf {input} --make-bed --out CRA_annotated_plink_chr{params.chrom}
		"""
		
rule merge:
	input:
		expand("CRA_annotated_plink_chr{chr}.fam",chr=CHRS),
		expand("CRA_annotated_plink_chr{chr}.bim",chr=CHRS),
		expand("CRA_annotated_plink_chr{chr}.bed",chr=CHRS)
	output:
		bed="CRA_annotated_plink_merged.bed",
		bim="CRA_annotated_plink_merged.bim",
		fam="CRA_annotated_plink_merged.fam"
	shell:
		"""
		mplink19 --bfile CRA_annotated_plink_chr1 --merge-list bmergelist --geno 0.02 --maf 0.0001 --make-bed --out CRA_annotated_plink_merged
		"""	

rule get_pedigree:
	input:
		"CRA_annotated_plink_merged.fam"
	output:
		"CRA_annotated_plink_merged_ped.fam"
	shell:
		"""
		Rscript put_in_pedigree.R
		"""
	
rule sex_check:
	input:
		bed="CRA_annotated_plink_merged.bed",
		bim="CRA_annotated_plink_merged.bim",
		fam="CRA_annotated_plink_merged.fam"
	output:
		"sexcheck.sexcheck"
	shell:
		"""
		mplink19 --bed {input.bed} --bim {input.bim} --fam {input.fam} --check-sex --out sexcheck
		"""		

rule king:
	input:
		bed="CRA_annotated_plink_merged.bed",
		bim="CRA_annotated_plink_merged.bim",
		fam="CRA_annotated_plink_merged.fam"
	output:
		"king_duplicate.con"
	shell:
		"""
		king -b {input.bed} --fam {input.fam} --bim {input.bim} --duplicate --prefix king_duplicate
		"""	

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.chromosomes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "XY"]
params.outputDir = "."
params.basePath = "/proj/regeps/regep00/studies/TopMed/data/dna/whole_genome/TopMed/data/freezes/freeze.10a"
params.samplesFile = "${params.basePath}/manifests/nwdids.txt"
params.minDP = "minDP10"

// These are given in GRCh38 coordinates
PAR1_RANGE = "chrX:10001-2781479"
PAR2_RANGE = "chrX:155701383-156030895"

process extract {
    conda "topmed_wgs_extraction/envs/bcftools.yaml"
    publishDir params.outputDir, mode: "copy"

    input:
    val chrom
    path bcf
    path samplesFile
    
    output:
    path "freeze.10a.chr${chrom}.${params.minDP}.bcf"

    
    script:
    output = "freeze.10a.chr${chrom}.${params.minDP}.bcf"
    if (chrom == 'X') {
        """
        bcftools view $bcf \
                --targets ^$PAR1_RANGE,$PAR2_RANGE \
                -s $samplesFile \
                -i 'FILTER=\"PASS\"' \
                -c 1 \
                --force-samples \
            | python ../scripts/python/chrX_remove_het_male_snps.py $samples 
            | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -O b -o $output
        """
	}
    else if (chrom == 'XY') {
        """
        bcftools view $bcf \
                --regions $PAR1_RANGE,$PAR2_RANGE \
                -s $samplesFile \
                -i 'FILTER=\"PASS\"' \
                --force-samples \
                -c 1 \
            | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -O b -o $output
        """
	}
    else {
        """
        bcftools view -s $samplesFile -i 'FILTER=\"PASS\"' -c 1 --force-samples $bcf \
            | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -O b -o $output
        """
	}
}

process index {
    conda "topmed_wgs_extraction/envs/bcftools.yaml"
    publishDir params.outputDir, mode: "copy"

    input:
    path bcf
    
    output:
    path "${bcf}.csi"
    
    """
    bcftools index $bcf --force --csi
    """
}

workflow {

    // Prepare the input channels to match the run params
    input_bcf_files = []
    params.chromosomes.each {
        input_bcf_files = input_bcf_files + "${params.basePath}/${params.minDP}/freeze.10a.chr${it}.pass_and_fail.gtonly.${params.minDP}.bcf"
    }
    chromosomes = Channel.fromList(params.chromosomes)
    input_bcf_files = Channel.fromList(input_bcf_files)
    samples_file = Channel.of(params.samplesFile)

    // Run the workflow
    extract(chromosomes, input_bcf_files, samples_file)
    index(extract.output)
}

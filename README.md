# TOPMed WGS Freeze.10 Extraction

## ⚡Quickstart (In-Development)
```
$ /bin/bash
$ source /proj/relibs/relib00/conda-cdnm/bin/activate /udd/rejpz/.conda/envs/nextflow
$ topmed-wgs-extract --version
$ topmed-wgs-extract --pepfile=pep.yaml \
    --extract-dir=/d/tmp/regeps/regep00/studies/COPDGene/analyses/rejpz/
$ topmed-wgs-extract --configfile=config.yaml
```

## Overview
For several important reasons (e.g. need to combine datasets for downstream analysis, keep QC criteria the same across TOPMed, keep all variants for rare variant analysis), we have decided to keep the raw genotype data in the original format delivered by TopMed (QC performed by IRC); and require that (groups of) downstream analysts extract a working subset of samples and positions that are relevant to their specific analyses.

To this end, the Freeze.10 datasets are being released as a series of PEP datasets and a reproducibility tool that can (optionally) be used to extract the PEP dataset into a working location. [[ The PEP Datasets, Dataset Metadata are contained in a project-level yaml file, and the Sample Metadata are contained in the PEP Sample Table. Sample Metadata include a list of samples that have passed internal QC for a particular study (QC performed by CDNM), and associated metadata from both the Study (Subject Id, Aliases) and the Laboratory (Sample Type, Collection Information). ]]

## Configuration Options
Which samples to extract?
‼️ pepfile

The path to a valid, relevant PEP file (e.g. one from the table above under Canonical Dataset Locations). By default, the tool with use the “action” column of the SampleTable to extract only those samples with a value of “keep” in this column (override?)

Which directory to extract to?
‼️ extract_dir

You are encouraged to extract into the scratch space if you are familiar with its use (and if there is space). Please talk to Bioinformatics if you have not used the scratch space before. Exceptions can be made for data that must sit on the filesystem for an extended period of time.

Which files to extract from?
raw_wgs_base_path - pep.config['base_path']
depth - minDP0 minDP10

Which loci to extract?
```
bcftools view {input.bcf} \
    -s {params.samples} \
    -i 'FILTER="PASS"' \
    -c 1 \
    -O b \
    -m2 \
    -M2 \
    --force-samples \
    --types snps \
    -o {output.bcf}
```

## Output

## Discussion
This is certainly a departure from our normal data release process.

The Metadata in the PEP Datasets below can be used to drive this extraction: the paths to files from which data should be extracted; and a table of relevant passing samples to extract. Additional filters can be added as needed for variant-type, etc. These analytical datasets should be created in the scratch-space, and removed when they are no longer needed, so it is very important that the code used to create these datasets is recorded. Larger extractions can be slow, so make sure you put some thought into what is needed downstream. Extraction will benefit greatly from parallelization, which can take many forms depending on your level of comfort with the cluster.

http://pep.databio.org/en/2.1.0/specification/
http://pep.databio.org/


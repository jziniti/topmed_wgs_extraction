#!/bin/bash
source /proj/relibs/relib00/conda-cdnm/bin/activate sm6
snakemake -s /udd/rejpz/code/topmed_wgs_extraction/workflows/pep_extraction_example.smk \
          --config pepfile="LTRC_freeze.10_20220621.pep.yaml" \
          --cluster "qsub -v PATH -cwd -terse -S /bin/bash" \
          --latency-wait 120 \
          --use-conda \
          --conda-prefix "/proj/relibs/relib00/smk-conda-cache/envs/" \
          --jobs=23 \
	  -n


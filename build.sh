mkdir ${PREFIX}/snakemake/

git clone git@changit.bwh.harvard.edu:rejpz/topmed_wgs_extraction.git
cd topmed_wgs_extraction
cp -R envs workflows scripts notebooks ${PREFIX}/snakemake/

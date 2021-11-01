package_name=topmed_wgs_extraction
mkdir ${PREFIX}/snakemake/
mkdir ${PREFIX}/snakemake/${package_name}

git clone git@changit.bwh.harvard.edu:rejpz/${package_name}.git
cd ${package_name}
for dir in envs workflows scripts notebooks rules; do
    if [ -e ${dir} ]; then
        cp -R ${dir} ${PREFIX}/snakemake/${package_name}/
    fi
done

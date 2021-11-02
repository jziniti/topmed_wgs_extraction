import pandas as pd
from pathlib import Path

def set_basename(row):
    abspath = Path(row['abspath'])
    return abspath.name

if __name__ == '__main__':
    bams_with_md5s = pd.read_csv(snakemake.input[0], delim_whitespace=True, names=('md5sum', 'abspath'))
    bams_with_md5s['basename'] = bams_with_md5s.apply(set_basename, axis=1)
    bams_with_md5s.to_csv(snakemake.output[0], index=False, header=False, columns=('abspath', 'basename'), sep="\t")

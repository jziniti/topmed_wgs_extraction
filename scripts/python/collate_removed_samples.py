import pandas as pd

if __name__ == '__main__':
    manifest = pd.read_csv(snakemake.input.manifest)
    dropped = manifest[(manifest['dataset_id'] == 0) & (manifest['action'] == 'drop')][['SampleAlias']]
    print(dropped)

    qc_flags = pd.read_csv(snakemake.input.flags, names=('_index', 'SampleAlias', 'reason'))
    qc_flags = qc_flags[['SampleAlias']]
    print(qc_flags)

    dropped_samples = pd.concat([dropped, qc_flags])
    dropped_samples['FID'] = '0'
    print(dropped_samples)
    dropped_samples.to_csv(snakemake.output.csv)
    dropped_samples.to_csv(snakemake.output.txt, sep='\t', index=False, header=False, columns=['FID', 'SampleAlias'])
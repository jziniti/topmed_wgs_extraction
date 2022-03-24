import pandas as pd
import logging

if __name__ == '__main__':

    log = logging.getLogger('topmed_wgs_extraction.collate_remove_samples')
    logging.basicConfig(level=logging.DEBUG)

    manifest = pd.read_csv(snakemake.input.manifest)
    dropped = manifest[(manifest['dataset_id'] == 0) & (manifest['action'] == 'drop')][['SampleAlias', 'flags']]
    dropped.rename(columns={'flags':'reason'}, inplace=True)
    log.debug(f'{dropped=}')

    automated_flags = pd.read_csv(snakemake.input.flags, index_col=0)
    automated_flags.rename(columns={'IID':'SampleAlias', 'Reason':'reason'}, inplace=True)
    log.debug(f'{automated_flags=}')

    dropped_samples = pd.concat([dropped, automated_flags])
    dropped_samples = dropped_samples.groupby('SampleAlias').agg({'reason': ';'.join}).reset_index()

    dropped_samples['FID'] = '0'
    dropped_samples.rename(columns={'SampleAlias':'sample_name'}, inplace=True)
    log.debug(f'{dropped_samples=}')
    dropped_samples.to_csv(snakemake.output.csv, index=False, columns=['sample_name', 'reason'])
    dropped_samples.to_csv(snakemake.output.txt, sep='\t', index=False, header=False, columns=['FID', 'sample_name']) ### This one is used by plink to do the removal
import pandas as pd

if __name__ == '__main__':

    manifest = pd.read_csv(snakemake.input.manifest)
    manifest['S_SAMPLEID'] = manifest.apply(lambda x: x['SampleLabel'].split('_')[0], axis=1)
    manifest.to_csv('tmp/resta_manifest.csv')

    king_results = []
    for file in snakemake.input.king_results:
        df = pd.read_csv(file, names=('FID1','ID1','FID2','ID2','N_SNP','HetHet','IBS0','Kinship'))
        df = df[df['FID1'] != 'FID1']
        king_results.append(df)
    all_king_results = pd.concat(king_results)
    all_king_results = all_king_results.merge(manifest, left_on='ID1', right_on='S_SAMPLEID', how='left')
    all_king_results.to_csv(snakemake.output.df)
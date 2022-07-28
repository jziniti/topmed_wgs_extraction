import pandas as pd

if __name__ == '__main__':

    masterfile = pd.read_csv(snakemake.input.masterfile, sep='\t')
    batch_set = masterfile[masterfile['batch'] == snakemake.wildcards.batch]
    if batch_set.shape[0] == 0:
        print(f'No sample from the masterfile match batch {snakemake.wildcards.batch}, using whole mastefile')
        batch_set = masterfile

    master_batch = batch_set[batch_set['final.analyzed'] == 1]
    print(master_batch)
    
    ### The batch data needs to be modified to conform to the masterfile
    batch = pd.read_csv(snakemake.input.batch, usecols=['IID_orig_x', 'IID_orig_y', 'Kinship'])  #FIXME: This needs to be changed upstream
    batch_rna_wgs = batch[(batch.IID_orig_x.str.startswith('NWD') & batch.IID_orig_y.str.startswith('S-'))]
    batch_rna_wgs['sample_id'] = batch_rna_wgs.apply(lambda row: row['IID_orig_y'].split('_')[0], axis=1)
    print(batch_rna_wgs)
    batch_rna_wgs = batch_rna_wgs.merge(master_batch, on='sample_id', how='left')
    
    print(batch_rna_wgs)

    groups_by_sample = batch_rna_wgs.groupby('BAM_file')
    filtered_results = []
    print(groups_by_sample.groups)
    for (rna_sample, group) in groups_by_sample:
        max = group['Kinship'].max()
        s25 = max * 0.25
        s50 = max * 0.50
        s90 = max * 0.90
        n0 = group[group.Kinship < s25].size
        n25 = group[(group.Kinship >= s25) & (group.Kinship < s50)].size
        n50 = group[(group.Kinship >= s50) & (group.Kinship < s90)].size
        tranch90 = group[group.Kinship >= s90]
        n90 = tranch90.shape[0]
        unique_hits = list(set(tranch90['IID_orig_x'].values))
        if len(unique_hits) <= 2:
            filtered_results.append(tranch90)
        print(f'{rna_sample} {max}: {s25}/{s50}/{s90} {n0}-{n25}-{n50}-{n90} ({tranch90})')
    
        ### Save the n75 group to the filtered_results
    if master_batch.shape[0] == 0:
        filtered_results = pd.DataFrame(columns=batch_rna_wgs.columns)
    else:
        print(filtered_results)
        filtered_results = pd.concat(filtered_results)
    filtered_results.to_csv(snakemake.output.csv, index=False)

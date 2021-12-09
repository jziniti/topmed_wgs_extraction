import pandas as pd

if __name__ == '__main__':

    ### Read teh pheno file to get gender
    pheno = pd.read_csv(snakemake.input.phenofile, sep='\t')

    ### Read and process the input
    masterfile = pd.read_csv(snakemake.input.masterfile, sep='\t')
    final_analyzed = masterfile[masterfile['final.analyzed'] == 1]
    final_analyzed = final_analyzed.merge(pheno, left_on='COPDgeneID', right_on='sid', how='left')
    final_analyzed['ObservedGender'] = final_analyzed['gender'].replace({1.0: 'M', 2.0: 'F'})
    rename_columns = {
        'BAM_file':'SampleLabel',
        #'sample':'SampleAlias',
        'actual_id':'SampleAlias',
        'sample_id':'S_SAMPLEID',
        }
    final_analyzed.rename(columns=rename_columns, inplace=True)

    ### Write a filtered output manifest
    final_analyzed.to_csv(snakemake.output.manifest, index=False, columns=['SampleLabel', 'SampleAlias', 'batch',]) #'actual_id',#'S_SAMPLEID'])

    ### Gender Observations
    final_analyzed.to_csv(snakemake.output.gobs, index=False, columns=['SampleLabel', 'ObservedGender',])

    ### Rassignments
    ## reassignments = final_analyzed[final_analyzed['COPDgeneID'] != final_analyzed['actual_id']]
    reassignments = final_analyzed[final_analyzed['COPDgeneID'] != final_analyzed['SampleAlias']]
    reassignments.to_csv(snakemake.output.reassignments, columns=['SampleLabel', 'SampleAlias'])
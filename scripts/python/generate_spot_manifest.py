### Generate something that can be loaded into spot
### S_SAMPLEID, SDI_ALIASID, ITEMID, PATH, QC_STATUS, QC_DESCRIPTION

import pandas as pd

def set_qc_status(row):
    keep_drop = row['keep_drop']
    if  keep_drop == 'keep':
        qc_status = 'Pass'
    else:
        qc_status = 'InternalQcFail'
    return qc_status

if __name__ == '__main__':

    s_studyid = snakemake.wildcards.s_studyid
    study_metadata = snakemake.config[s_studyid]
    dataset_ids = []
    dataset_paths = []
    load_into_spot = []
    for dataset in study_metadata['datasets']:
        id = dataset['id']
        manifest_path = dataset['manifest_path']
        path = dataset.get('spot_path', manifest_path)
        load = dataset.get('load_into_spot', False)
        dataset_ids.append(id)
        dataset_paths.append(path)
        load_into_spot.append(load)
    dataset_metadata = pd.DataFrame({'dataset_id':dataset_ids, 'PATH':dataset_paths, 'load_into_spot': load_into_spot})
    print(dataset_metadata)

    multiomics_manifest = pd.read_csv(snakemake.input.manifest)
    multiomics_manifest = multiomics_manifest.merge(dataset_metadata, on='dataset_id', how='left')
    multiomics_manifest = multiomics_manifest[multiomics_manifest['load_into_spot'] == True]
    print(multiomics_manifest)

    multiomics_manifest['SDI_ALIASID'] = multiomics_manifest['SampleAlias']
    multiomics_manifest['ITEMID'] = multiomics_manifest['SampleLabel']
    try:
        multiomics_manifest['QC_STATUS'] = multiomics_manifest.apply(set_qc_status, axis=1)
    except:
        multiomics_manifest['QC_STATUS'] = None
    multiomics_manifest['QC_DESCRIPTION'] = multiomics_manifest['status']
    print(multiomics_manifest)
    multiomics_manifest.to_csv(snakemake.output.manifest, index=False, columns=('S_SAMPLEID', 'SDI_ALIASID', 'ITEMID', 'PATH', 'QC_STATUS', 'QC_DESCRIPTION'))
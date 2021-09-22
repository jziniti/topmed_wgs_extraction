import pandas as pd

if __name__ == '__main__':
    manifest = pd.read_csv(snakemake.input.manifest)
    manifest_columns = manifest.columns

    plate_manifests = []
    for filename in snakemake.input.plate_manifests:
        plate_manifest = pd.read_csv(filename)
        plate_manifests.append(plate_manifest)
    all_plate_103 = pd.concat(plate_manifests)

    plate_103_subset = manifest.merge(all_plate_103, how='right', left_on='S_SUBJECTID', right_on='SUBJECTID')
    group_ids = plate_103_subset.group_id.unique()

    manifest = manifest[manifest['group_id'].isin(group_ids)]
    manifest.to_csv(snakemake.output.manifest)
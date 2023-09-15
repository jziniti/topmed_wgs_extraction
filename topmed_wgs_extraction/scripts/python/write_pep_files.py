"""
http://pep.databio.org/en/latest/
http://pep.databio.org/en/latest/specification/

"""

import getpass
import logging
import os
import pandas as pd
import pathlib
import yaml


SAMPLE_TABLE_COLUMNS = ['sample_name', 'S_SAMPLEID', 'S_SUBJECTID', 'S_STUDYID', 'ALIASTYPE', 'ALIASID', 'batchid',
                        's_samplefamilyid', 'tissueid', 'sampletypeid', 'u_initsampletypeid', 'u_collectionid', 'u_specimentypeid', 'ExpectedGender',
                        'ObservedGender', 'flags', 'action', 'reason']

def default_keep(row):
    override_actions = {'', 'reassign'}
    current_action = row['action']
    action = current_action
    if pd.isna(current_action) or current_action in override_actions:
        action = 'keep'
    return action

def default_drop(row):
    allowed_actions = {'keep',}
    current_action = row['action']
    action = current_action
    if pd.isna(current_action) or current_action not in allowed_actions:
        action = 'drop'
    return action

def label_duplicates(row):
    current_action = row['action']
    sample_name = row['SampleLabel']
    stid = row['S_SUBJECTID']
    all_samples_for_subject = sorted(list(manifest[(manifest['S_SUBJECTID'] == stid) & (manifest['action'] == current_action)]['SampleLabel']))
    if current_action == 'keep' and sample_name != all_samples_for_subject[0]:
        action = 'duplicate'
    else:
        action = current_action
    return action

if __name__ == '__main__':
    log = logging.getLogger('topmed_freeze10.write_pep_files')
    logging.basicConfig(level=logging.DEBUG)

    snakemake = snakemake # type: ignore
    with open(snakemake.output.pep, 'w') as pm_file:
        files = []
        wgs_path_wildcard = snakemake.params['raw_wgs_path_wildcard']
        (base_path, filename_wildcard) = wgs_path_wildcard.rsplit('/', 1)
        for chrom in list(snakemake.params.chromosomes):
            files.append(filename_wildcard.format(chrom=chrom))
        pep_project = {
            'pep_version': '2.0.0',
            'sample_table': str(pathlib.Path(snakemake.output.sample_table).name),
            'chromosomes': list(snakemake.params.chromosomes),
            'split_by_chromosome': True,
            # 'files': list([pathlib.Path(x).name for x in snakemake.input.bcfs]),
            'files': files,
            'cwd': os.getcwd(),
            'user': getpass.getuser(),
            'cmdline': snakemake.params.get('cmdline', ''),
            'timestamp': snakemake.params.get('timestamp', ''),
            'reference_genome': snakemake.params.get('reference_genome', ''),
            'issue_url': snakemake.params.get('issue_url', ''),
            'base_path': base_path,
        }
        pm_file.write(yaml.dump(pep_project))

    ### Read the metrics file (we will use this to drop/keep the correct sample)
    sample_annotations = pd.read_csv(snakemake.input.sample_annotations, delim_whitespace=True)
    log.debug(f'{sample_annotations=}')

    manifest = pd.read_csv(snakemake.input.manifest)
    manifest = manifest[manifest['dataset_id'] == 0] ## FIXME: this isn't always true ...
    manifest = manifest.merge(sample_annotations, left_on='SampleAlias', right_on='SAMPLE_ID', how='left')
    
    manifest['action'] = manifest.apply(default_keep, axis=1)
    
    ### Flag duplicates
    manifest = manifest.sort_values(['SampleAlias', 'FRAC_DP20'], ascending=False) ## Sort so that the highest FRAC_DP20 is kept
    manifest['action'] = manifest.apply(label_duplicates, axis=1) ## These are labeled, only to be dropped next step

    manifest = manifest[manifest['action'] == 'keep']
    manifest.rename(columns={'SampleLabel':'sample_name', 'sampletypeid_x': 'sampletypeid'}, inplace=True)
    manifest.to_csv(snakemake.output.sample_table, index=False, columns=SAMPLE_TABLE_COLUMNS)

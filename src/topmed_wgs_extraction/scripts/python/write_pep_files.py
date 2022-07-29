"""
http://pep.databio.org/en/latest/
http://pep.databio.org/en/latest/specification/

"""

import pandas as pd
import pathlib
import yaml

import os
import getpass

SAMPLE_TABLE_COLUMNS = ['sample_name', 'S_SAMPLEID', 'S_SUBJECTID', 'S_STUDYID', 'ALIASTYPE', 'ALIASID', 'batchid',
                        's_samplefamilyid', 'tissueid', 'sampletypeid_x', 'sampletypeid_y', 'u_collectionid', 'u_specimentypeid', 'ExpectedGender',
                        'ObservedGender', 'flags', 'group_id', 'c_group', 'action']

if __name__ == '__main__':
    with open(snakemake.output.pep, 'w') as pm_file:
        pep_project = {
            'pep_version': '2.0.0',
            'sample_table': str(pathlib.Path(snakemake.output.sample_table).name),
            'chromosomes': list(snakemake.params.chromosomes),
            'split_by_chromosome': True,
            'files': list([pathlib.Path(x).name for x in snakemake.input.bcfs]),
            'cwd': os.getcwd(),
            'user': getpass.getuser(),
            'cmdline': snakemake.params.get('cmdline', ''),
            'timestamp': snakemake.params.get('timestamp', ''),
            'reference_genome': snakemake.params.get('reference_genome', ''),
            'issue_url': snakemake.params.get('issue_url', ''),
        }
        pm_file.write(yaml.dump(pep_project))

    manifest = pd.read_csv(snakemake.input.manifest)
    manifest = manifest[manifest['dataset_id'] == 0] ## FIXME: this isn't alwayts true ...
    manifest.rename(columns={'SampleLabel':'sample_name'}, inplace=True)
    
    manifest.to_csv(snakemake.output.sample_table, index=False, columns=SAMPLE_TABLE_COLUMNS)

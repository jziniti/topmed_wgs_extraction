#!/usr/bin/env python3

import argparse
from datetime import date, datetime
import logging
import os
from pathlib import Path
import sys
from shutil import copyfile
import snakemake

#VERSION = '2.3'
PACKAGE_DIR = Path(os.path.dirname(__file__)).absolute()

with open(PACKAGE_DIR/"VERSION",'r') as versionFile:
    VERSION = versionFile.readline().strip()
    
PROFILE_NAME = 'cdnm'

def main():
    logging.basicConfig(level='INFO')
    LOG = logging.getLogger('topmed_wgs_extraction.run_workflow.main')
    LOG.info('topmed_wgs_extraction.run_workflow.main()')
    
    workflows_path = PACKAGE_DIR/'rules'
    profile_dir = PACKAGE_DIR/'profiles'/PROFILE_NAME 
    workflow_path = workflows_path/'pep_extraction_example.smk'

    #setup argparser to display help if no arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description=f"Extract TOPMed WGS Data for Downstream Analysis", usage=f"topmed-wgs-extract --configfile=<configfile>")
    parser.add_argument('-v', '--version', action='store_true', help="Print Version Number")
    parser.add_argument('--configfile', type=str, help="The configuration file containing the run settings")
    parser.add_argument('--pepfile', type=str, help="The PEP Dataset containing the samples to extract")
    parser.add_argument('--extract-dir', type=str, help="The directory location to extract the files into")
    parser.add_argument('--get-config', action="store_true", help="Get the configuration template for this workflow")

    args = parser.parse_args()

    config = {
        'pepfile':args.pepfile,
        'extract_dir': args.extract_dir,
        }
    configfiles = []
    if args.configfile:
        configfiles.append(args.configfile)
    
    #give config to user if requested
    if args.version:
        print(f'{VERSION}')
        sys.exit()
    elif args.get_config:
        cwd = Path(os.getcwd())
        config_path = PACKAGE_DIR/'conf/sample.config.dn8.yaml'
        #config_data = files('conf').joinpath('sample.config.dn8.yaml').read_text()
        datestamp = date.today().strftime('%Y-%m-%d')
        dest_path = cwd/f"config-civic-{datestamp}.yaml"
        LOG.info(f'{config_path=} => {dest_path=}')
        #with open(dest_path, 'w') as dst:
        #    dst.write(config_data)
        copyfile(config_path, dest_path)
        sys.exit()

    #LOG.info('Running the workflow')
    #workdir
    #run_id
    #config-timestamp
    #stats
    #summary
    run_id = datetime.now().strftime('%Y%m%d_%H%M%S')
    snakemake.snakemake(workflow_path,
                        configfiles=configfiles,
                        config=config,
                        use_conda=True,
                        printshellcmds=True,
                        cluster="qsub -v PATH -cwd -l lx -terse -S /bin/bash",
                        latency_wait=120,
                        conda_prefix="/proj/relibs/relib00/smk-conda-cache/envs/",
                        printreason=True,
                        nodes=25,
                        workdir=f'/d/tmp2/log/rejpz/multiomics_qc/{run_id}/',
    )

if __name__ == '__main__':
    main()
    

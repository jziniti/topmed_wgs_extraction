#!/usr/bin/env python

import os
import setuptools
import versioneer
 
PACKAGE_NAME = 'topmed_wgs_extraction'

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

package_data = package_files(f'{PACKAGE_NAME}/conf') + package_files(f'{PACKAGE_NAME}/rules') + \
               package_files(f'{PACKAGE_NAME}/envs') + package_files(f'{PACKAGE_NAME}/scripts')
package_data.append(f'../{PACKAGE_NAME}/VERSION')

# VERSION = '1.0'
#with open(f"{PACKAGE_NAME}/VERSION",'r') as vf:
#    (major_minor, point, hash) = vf.readline().strip().split('-')
#    VERSION = f'{major_minor}.{point}'

setuptools.setup(
    name=PACKAGE_NAME,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="John Ziniti",
    author_email="john.ziniti@channing.harvard.edu",
    description="Code/package for extracting data from TOPMed Freeze10 WGS",
    url="https://changit.bwh.harvard.edu/rejpz/topmed_wgs_extraction",
    long_description="",
    long_description_content_type="text/markdown",

    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={
        '': package_data
    },

    entry_points={
        "console_scripts": [
            'topmed-wgs-extract = topmed_wgs_extraction.run_workflow:main'
        ]
    },
)

#!/usr/bin/env python

import os
import setuptools

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

package_data = package_files('topmed_wgs_extraction/conf') + package_files('topmed_wgs_extraction/rules') + package_files('topmed_wgs_extraction/envs') + package_files('topmed_wgs_extraction/scripts')

VERSION = '1.0'

setuptools.setup(
    name="topmed_wgs_extraction",
    version=VERSION,
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

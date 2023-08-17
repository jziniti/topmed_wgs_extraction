#!/bin/env python
"""
Tried running the chrX data on UMich server and it failed with:


Chromosome X check failed! 
java.io.IOException: Found haplotype 0/1 at pos 2700157 for male proband 1674_1674
Found haplotype 0/1 at pos 2777560 for male proband 1434_1434
Found haplotype 0/1 at pos 2777560 for male proband 1937_1937
Found haplotype 0/1 at pos 2777560 for male proband 3711_3711
...

etc.  

Apparently, heterozygous SNPs for males need to be removed.
"""

import gzip
import logging
import pandas as pd
import sys

HAPLOID_CALLS = {
    './.': '.', 
    '0/1': '.', 
    '1/0': '.',
    '0/0': '0',
    '1/1': '1',
    }

if __name__ == '__main__':
    LOG = logging.getLogger(__name__)
    logging.basicConfig(level=logging.DEBUG)

    #gender_column = snakemake.params.get('gender_column', 'ExpectedGender')
    #male_gender_value = snakemake.params.get('male_gender_value', 'M')
    #sample_table = snakemake.params.sample_table
    gender_column = 'ExpectedGender'
    male_gender_value = 'M'
    sample_table = pd.read_csv(sys.argv[1])
    males = set(list(sample_table[sample_table[gender_column] == male_gender_value]['sample_name']))
    LOG.debug(f'{len(males)=}')

    dst = sys.stdout
    for line in sys.stdin:
        if line.startswith('##'):
            dst.write(line)
        elif line.startswith('#CHROM'):
            fields = line.strip().split()
            male_columns = []
            for (f, field) in enumerate(fields):
                if field in males:
                    male_columns.append(f)
            dst.write(line)
            LOG.debug(f'{len(male_columns)=}')
        else:
            fields = line.strip().split()
            LOG.debug(f'Processing {fields[0-2]=}')
            # fields[0] = 'X'
            for f in male_columns:
                field = fields[f]
                fields[f] = HAPLOID_CALLS[field]
            line = '\t'.join(fields)
            dst.write(line)
            dst.write('\n')
                

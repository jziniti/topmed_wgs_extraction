"""
http://pep.databio.org/en/latest/
http://pep.databio.org/en/latest/specification/

"""

## import peppy
import pandas

## subsample_table: ["path/to/subsample_table.csv", "path/to/subsample_table2.csv"]

if __name__ == '__main__':
    with open(snakemake.output.project_metadata, 'w') as pm_file:
        pm_file.write(f"""pep_version: 2.0.0
sample_table: "{snakemake.output.sample_table}"
sample_modifiers:
  derive:
    attributes: [read1, read2, other_attr]
    sources:
      key1: "path/to/derived/value/{attribute1}"
      key2: "path/to/derived/value/{attr2}"
"""
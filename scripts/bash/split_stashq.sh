#!/bin/bash
INPUT_FILE=${1}
PREFIX=`echo $RANDOM | md5sum | head -c 20`

split -l 1000 ${INPUT_FILE} tmp/${PREFIX}.
for split in `ls tmp/${PREFIX}.*`; do
    stashq alias `cat ${split}`
done

INPUT_PREFIX=${1/.bed/}
OUTPUT_PREFIX=${2/.bed/}
plink --bfile $INPUT_PREFIX --out $OUTPUT_PREFIX --extract $3 --update-chr $4 --make-bed 
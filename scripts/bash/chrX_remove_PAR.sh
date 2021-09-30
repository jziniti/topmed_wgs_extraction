INPUT_PREFIX=${1/.bed/}
OUTPUT_PREFIX=${2/.bed/}
plink --bfile $INPUT_PREFIX --out $OUTPUT_PREFIX --exclude $3 --make-bed 
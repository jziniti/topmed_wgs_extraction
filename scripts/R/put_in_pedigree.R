library(tidyverse)

fam=read.table(snakemake@input[["fam"]],as.is=T)
## "/proj/regeps/regep00/studies/CRA/metadata/CRA.fam"
ped=read.table(snakemake@input[["canonical_fam"]],as.is=T)

#"/proj/regeps/regep00/studies/CRA/metadata/cra_nwd_map.txt"
translation=read.table(snakemake@input[["sample_mapping"]],as.is=T)

print(dim(translation))
print(length(unique(translation[,1])))
print(length(unique(translation[,2])))

tmp=left_join(fam, translation, by=c("V2"="V2"))
tmp2=left_join(tmp, ped, by=c("V1.y"="V2"))


fam_out=tmp2[,c(8,7,9,10,11,12)]

write.table(fam_out, file=snakemake@output[["fam"]], col.names=F, row.names=F, quote=F)

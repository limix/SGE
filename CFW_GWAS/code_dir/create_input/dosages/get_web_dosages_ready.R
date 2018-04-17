chr <- as.numeric(commandArgs(trailingOnly = T)[1])

load('data_dir/dosages/web_dosages/chr19.prunedgen.final.maf001.0.98.RData')
all_mice = unlist(nameList)
rm(list = c('pruned_pos','pruned_dosages','nameList'))

print(chr)
load(paste('data_dir/dosages/web_dosages/chr',chr,'.prunedgen.final.maf001.0.98.RData',sep=''))
pruned_pos$ID = paste(as.character(pruned_pos[,'CHR']),pruned_pos[,'POS'],sep='_')
dosages = t(pruned_dosages)
colnames(dosages) = pruned_pos$ID
if (any(unlist(nameList)!=all_mice)) stop('pn mouse IDs')
all_mice = sub('_recal.reheadered.bam','',all_mice,fixed = T)
all_mice = sub('___','/',all_mice)
rownames(dosages) = all_mice

save(dosages,file=paste('data_dir/dosages/pruned_dosages/unstd_my_final_dosages_chr',chr,'.RData',sep=''))

std_dosages = scale(dosages,center = T, scale = T)
save(std_dosages,file=paste('data_dir/dosages/pruned_dosages/std_my_final_dosages_chr',chr,'.RData',sep=''))

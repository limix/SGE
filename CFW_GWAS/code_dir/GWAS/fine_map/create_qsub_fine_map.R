load('output_dir/pvalues_pointwise_conditional/SGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData')
all_QTLs = all_QTLs[order(all_QTLs[,'chr'],all_QTLs[,'pos']),]

toWrite = NULL
#for (i in 1:dim(all_QTLs)[1]) {
for (i in c(2,3,7,13,16)) {

	toWrite = rbind(toWrite,paste('bash code_dir/GWAS/fine_map/fine_map.sh',all_QTLs[i,'measure'],all_QTLs[i,'chr'],all_QTLs[i,'pos'],all_QTLs[i,'ci_starts'],all_QTLs[i,'ci_stops']))
}

write.table(toWrite,file = 'code_dir/GWAS/fine_map/paper_qsub_fine_map.sh',quote = F, col.names = F, row.names = F)






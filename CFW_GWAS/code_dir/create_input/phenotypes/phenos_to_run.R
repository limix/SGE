load('data_dir/phenotypes/all_pheno_data.RData')
g = grep('_BWcorr$',colnames(phenotypes_bc))
remove  = sub('_BWcorr','',colnames(phenotypes_bc)[g])
w_remove = match(remove,colnames(phenotypes_bc))
keep = seq(1,dim(phenotypes_bc)[2])[-w_remove]
length(keep)
#[1] 170

write.table(colnames(phenotypes_bc),file = 'data_dir/phenotypes/colnames_phenotypes_bc.txt',col.names = F,row.names = F, sep='\t',quote = F)
colnames(phenotypes_bc)[keep]
write.table(colnames(phenotypes_bc)[keep],file = 'data_dir/phenotypes/colnames_phenotypes_bc_to_run.txt',col.names = F,row.names = F, sep='\t',quote = F)
write.table(keep,file = 'data_dir/phenotypes/phenos_to_run.txt',col.names = F,row.names = F, sep='\t',quote = F)

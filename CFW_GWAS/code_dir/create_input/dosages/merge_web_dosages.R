my_f = function(chr) {
	print(chr)
	load(paste('data_dir/dosages/pruned_dosages/unstd_my_final_dosages_chr',chr,'.RData',sep=''))
	load(paste('data_dir/dosages/pruned_dosages/std_my_final_dosages_chr',chr,'.RData',sep=''))
	return(list(dosages,std_dosages))
}
res = lapply(1:19,my_f)
dosages = do.call('cbind',lapply(res,'[[',1))
dim(dosages)
std_dosages = do.call('cbind',lapply(res,'[[',2))
dim(dosages)

load('data_dir/subsets/globally_excluded.RData')
motch = match(globally_excluded,rownames(dosages))
dosages = dosages[-na.exclude(motch),]
dim(dosages)
motch = match(globally_excluded,rownames(std_dosages))
std_dosages = std_dosages[-na.exclude(motch),]
dim(dosages)

save(dosages,file = 'data_dir/dosages/pruned_dosages/unstd_my_final_dosages.RData')
save(std_dosages,file = 'data_dir/dosages/pruned_dosages/std_my_final_dosages.RData')


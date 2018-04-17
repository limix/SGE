load('data_dir/dosages/pruned_dosages/GRM.RData')
load('data_dir/phenotypes/imputed_cage_info.RData')
load('data_dir/dosages/close_pairs_GRM.RData')

inter=intersect(colnames(GRM),rownames(imputed_cage_info))
inter = inter[!inter %in% all_cagemates_close_pairs]
length(inter)
#[1] 1812
GRM=GRM[match(inter,rownames(GRM)),match(inter,colnames(GRM))]
imputed_cage_info=imputed_cage_info[match(inter,rownames(imputed_cage_info)),]
cage_cm = imputed_cage_info[,'cage_WH']
all_mice = rownames(GRM)


all_mice = sample(all_mice,size = 1800, replace = F)
length(all_mice)
#[1] 1800 so that can be divided by 2,3,4,5 and 6

fake_cage_cm = matrix(NA, ncol = 5, nrow = length(all_mice))
colnames(fake_cage_cm) = paste('cage_gs',2:6,sep='')
rownames(fake_cage_cm) = all_mice
for (gs in 2:6) {
	nb_cages = length(all_mice)/gs
	cages = paste('cage',1:nb_cages,sep='')
	fake_cage_cm[,paste('cage_gs',gs,sep='')] = rep(cages, each = gs)
}

save(fake_cage_cm, file = 'data_dir/simulations/power_analysis/fake_cage_cm.RData')

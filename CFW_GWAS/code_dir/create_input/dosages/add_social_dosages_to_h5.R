load('data_dir/dosages/pruned_dosages/unstd_my_final_dosages.RData')

load('data_dir/phenotypes/imputed_cage_info.RData')
motch=match(rownames(dosages),rownames(imputed_cage_info))
any(is.na(motch))
#FALSE
imputed_cage_info = imputed_cage_info[motch,]

social_dosage_sumUnstd = matrix (NA,ncol = dim(dosages)[2],nrow = dim(dosages)[1]) 
rownames(social_dosage_sumUnstd) = rownames(dosages)
colnames(social_dosage_sumUnstd) = colnames(dosages)

for (mouse in rownames(dosages)) {
	first_cage = imputed_cage_info[mouse,1]
	cagemates = rownames(imputed_cage_info)[which(imputed_cage_info[,1]==first_cage)]
	cagemates = cagemates[cagemates!=mouse]

	if (length(cagemates) != 2) stop('pb')

	#no NA in dosages
	social_dosage_sumUnstd[mouse,] = colSums(dosages[cagemates,,drop=F])
}

any(is.na(social_dosage_sumUnstd))
#[1] FALSE

dim(social_dosage_sumUnstd)
#[1]   1869 353697

save(social_dosage_sumUnstd,file='data_dir/dosages/pruned_dosages/unstd_social_dosages_sumUnstd.RData')

std_social_dosage_sumUnstd = scale(social_dosage_sumUnstd,center = T, scale = T)
save(std_social_dosage_sumUnstd,file= 'data_dir/dosages/pruned_dosages/std_social_dosages_sumUnstd.RData')


load('data_dir/dosages/pruned_dosages/std_social_dosages_sumUnstd.RData')
library(rhdf5)
fid='data_dir/CFWmice.h5'

h5createGroup(fid,"social_pruned_dosages")

h5createGroup(fid,"social_pruned_dosages/row_header")
max(nchar(rownames(std_social_dosage_sumUnstd)))
#15
h5createDataset(file=fid,dataset="social_pruned_dosages/row_header/sample_ID",dims=dim(std_social_dosage_sumUnstd)[1],storage.mode='character',size=20)
h5write(obj=rownames(std_social_dosage_sumUnstd),file=fid,name="social_pruned_dosages/row_header/sample_ID")

pos = do.call('rbind',sapply(colnames(std_social_dosage_sumUnstd),strsplit,'_'))
pos = as.data.frame(pos,stringsAsFactors=F)
pos[,1] = as.numeric(sub('chr','',pos[,1]))
pos[,2] = as.numeric(pos[,2])
pos = as.matrix(pos)
ordre = order(pos[,1],pos[,2])
all(ordre == seq(1,dim(pos)[1]))
#TRUE
cumpos = pos[pos[,1]==1,2]
for (chr in unique(pos[,1])[-1]) {
	w = which(pos[,1]==chr)
	cumpos = c(cumpos,pos[w,2]+cumpos[length(cumpos)])
}
#plot(pos[,1],cumpos)
ordre = order(cumpos)
all(ordre == seq(1,length(cumpos)))
#TRUE 

h5createGroup(fid,"social_pruned_dosages/col_header")
h5write(obj=pos[,1],file=fid,name="social_pruned_dosages/col_header/chr")
h5write(obj=pos[,2],file=fid,name="social_pruned_dosages/col_header/pos")
h5write(obj=cumpos,file=fid,name="social_pruned_dosages/col_header/cumpos")
h5write(obj=colnames(std_social_dosage_sumUnstd),file=fid,name="social_pruned_dosages/col_header/id")

h5createDataset(file=fid,dataset="social_pruned_dosages/matrix",dims=dim(std_social_dosage_sumUnstd),chunk = c(dim(std_social_dosage_sumUnstd)[1],ceiling(dim(std_social_dosage_sumUnstd)[2])/2))
h5write(obj=std_social_dosage_sumUnstd,file=fid,name="social_pruned_dosages/matrix")



args <- commandArgs(TRUE)
locus = args[1]

load(paste('output_dir/fine_map_files/all_local_dosages/',locus,'.RData',sep=''))
dosages = t(dosages)
colnames(dosages) = map[,'id']
rownames(dosages) = sampleNames

load('data_dir/phenotypes/imputed_cage_info.RData')
motch=match(rownames(dosages),rownames(imputed_cage_info))
if (any(is.na(motch))) stop('pb')
imputed_cage_info = imputed_cage_info[motch,]

social_dosage_sumUnstd = matrix (NA,ncol = dim(dosages)[2],nrow = dim(dosages)[1]) 
rownames(social_dosage_sumUnstd) = rownames(dosages)
colnames(social_dosage_sumUnstd) = colnames(dosages)

for (mouse in rownames(dosages)) {
	first_cage = imputed_cage_info[mouse,1]
	cagemates = rownames(imputed_cage_info)[which(imputed_cage_info[,1]==first_cage)]
	cagemates = cagemates[cagemates!=mouse]
	social_dosage_sumUnstd[mouse,] = colSums(dosages[cagemates,,drop=F])
}

std_social_dosage_sumUnstd = scale(social_dosage_sumUnstd,center = T, scale = T)
dir.create('output_dir/fine_map_files/local_social_dosages')
save(std_social_dosage_sumUnstd,file = paste('output_dir/fine_map_files/local_social_dosages/',locus,'.RData',sep=''))

library(rhdf5)
fid='data_dir/CFWmice_locals.h5'
h5createFile(fid)

h5createGroup(fid,"local_social_dosages")
h5createGroup(fid,"local_social_dosages/row_header")
max(nchar(rownames(std_social_dosage_sumUnstd)))
#15
h5createDataset(file=fid,dataset="local_social_dosages/row_header/sample_ID",dims=dim(std_social_dosage_sumUnstd)[1],storage.mode='character',size=20)
h5write(obj=rownames(std_social_dosage_sumUnstd),file=fid,name="local_social_dosages/row_header/sample_ID")

h5createGroup(fid,paste("local_social_dosages/",locus,sep=""))
h5createGroup(fid,paste("local_social_dosages/",locus,"/col_header",sep=""))
h5write(obj=map[,'id'],file=fid,name=paste("local_social_dosages/",locus,"/col_header/id",sep=""))
h5write(obj=map[,'chr'],file=fid,name=paste("local_social_dosages/",locus,"/col_header/chr",sep=""))
h5write(obj=map[,'pos'],file=fid,name=paste("local_social_dosages/",locus,"/col_header/pos",sep=""))
h5write(obj=hwe_fail,file=fid,name=paste("local_social_dosages/",locus,"/col_header/hwe_fail",sep=""))
h5write(obj=info_fail,file=fid,name=paste("local_social_dosages/",locus,"/col_header/info_fail",sep=""))
h5write(obj=pruned_in,file=fid,name=paste("local_social_dosages/",locus,"/col_header/pruned_in",sep=""))

h5createDataset(file=fid,dataset=paste("local_social_dosages/",locus,"/matrix",sep=""),dims=dim(std_social_dosage_sumUnstd),chunk = c(dim(std_social_dosage_sumUnstd)[1],ceiling(dim(std_social_dosage_sumUnstd)[2])/2))
h5write(obj=std_social_dosage_sumUnstd,file=fid,name=paste("local_social_dosages/",locus,"/matrix",sep=""))

h5createGroup(fid,"local_direct_dosages/")
h5createGroup(fid,paste("local_direct_dosages/",locus,sep=""))
h5createDataset(file=fid,dataset=paste("local_direct_dosages/",locus,"/matrix",sep=""),dims=dim(dosages),chunk = c(dim(dosages)[1],ceiling(dim(dosages)[2])/2))
h5write(obj=dosages,file=fid,name=paste("local_direct_dosages/",locus,"/matrix",sep=""))

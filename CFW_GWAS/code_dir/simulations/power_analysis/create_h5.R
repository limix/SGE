load('data_dir/simulations/power_analysis/sims_design.RData')
	my_load = function(arg) {
    load(paste('data_dir/simulations/power_analysis/simulations/simulations_sd',arg,'.RData',sep=''))
    return(simulations)
}

res = lapply((dim(sims_design)[1] + 1) : (dim(sims_design)[1] + 24),FUN = my_load)

all_simulations = do.call('cbind',res)

all_simulations = all_simulations[,sample(1:dim(all_simulations)[2],replace = F)]
dim(all_simulations)

library(rhdf5)
fid='data_dir/simulations/power_analysis/simulations/simulations.hdf5'
h5createFile(fid)

h5createGroup(fid,"simulations")

h5createDataset(file=fid,dataset="simulations/array",dims=dim(all_simulations))
h5write(obj=all_simulations,file=fid,name="simulations/array")

h5createGroup(fid,"simulations/rows_subjects")
max(nchar(rownames(all_simulations)))
#15
h5createDataset(file=fid,dataset="simulations/rows_subjects/outbred",dims=dim(all_simulations)[1],storage.mode='character',size=20)
h5write(obj=rownames(all_simulations),file=fid,name="simulations/rows_subjects/outbred")

h5createGroup(fid,"simulations/cols_measures")
max(nchar(colnames(all_simulations)))
#31
h5createDataset(file=fid,dataset="simulations/cols_measures/measures",dims=dim(all_simulations)[2],storage.mode='character',size=35)
h5write(obj=colnames(all_simulations),file=fid,name="simulations/cols_measures/measures")


#dosages now. add genotypes for all SGE models in HDF5 file.
load('data_dir/simulations/power_analysis/dosages/dir_soc_std_dosages_gs3_colMax.RData')
std_dosages_gs3_colMax = std_dosages_gs3
std_social_dosage_sumUnstd_gs3_colMax = std_social_dosage_sumUnstd_gs3
for (gs in 2:6) {
	load(paste('data_dir/simulations/power_analysis/dosages/dir_soc_std_dosages_gs',gs,'.RData', sep=''))
}

std_dosages = cbind(std_dosages_gs2,std_dosages_gs3,std_dosages_gs3_colMax,std_dosages_gs4,std_dosages_gs5,std_dosages_gs6)
all_SNPs_to_save = unique(c(colnames(std_dosages_gs2),colnames(std_dosages_gs3),colnames(std_dosages_gs3_colMax),colnames(std_dosages_gs4),colnames(std_dosages_gs5),colnames(std_dosages_gs6)))
motch = match(all_SNPs_to_save, colnames(std_dosages))
any(is.na(motch))
#FALSE
std_dosages = std_dosages[,motch]
dim(std_dosages)
#[1]  1800 67,329

std_social_dosage_sumUnstd = cbind(std_social_dosage_sumUnstd_gs2,std_social_dosage_sumUnstd_gs3,std_social_dosage_sumUnstd_gs3_colMax,std_social_dosage_sumUnstd_gs4,std_social_dosage_sumUnstd_gs5,std_social_dosage_sumUnstd_gs6)
dim(std_social_dosage_sumUnstd)
#[1]  1800 10,1971

h5createGroup(fid,"direct_pruned_dosages")
h5createGroup(fid,"direct_pruned_dosages/row_header")
max(nchar(rownames(std_dosages)))
#15
h5createDataset(file=fid,dataset="direct_pruned_dosages/row_header/sample_ID",dims=dim(std_dosages)[1],storage.mode='character',size=20)
h5write(obj=rownames(std_dosages),file=fid,name="direct_pruned_dosages/row_header/sample_ID")
h5createGroup(fid,"direct_pruned_dosages/col_header")
max(nchar(colnames(std_dosages)))
#15
h5createDataset(file=fid,dataset="direct_pruned_dosages/col_header/id",dims=dim(std_dosages)[2],storage.mode='character',size=20)
h5write(obj=colnames(std_dosages),file=fid,name="direct_pruned_dosages/col_header/id")
h5createDataset(file=fid,dataset="direct_pruned_dosages/matrix",dims=dim(std_dosages),chunk = c(dim(std_dosages)[1],ceiling(dim(std_dosages)[2])/2))
h5write(obj=std_dosages,file=fid,name="direct_pruned_dosages/matrix")

h5createGroup(fid,"social_pruned_dosages")
h5createGroup(fid,"social_pruned_dosages/row_header")
max(nchar(rownames(std_social_dosage_sumUnstd)))
#15
h5createDataset(file=fid,dataset="social_pruned_dosages/row_header/sample_ID",dims=dim(std_social_dosage_sumUnstd)[1],storage.mode='character',size=20)
h5write(obj=rownames(std_social_dosage_sumUnstd),file=fid,name="social_pruned_dosages/row_header/sample_ID")
h5createGroup(fid,"social_pruned_dosages/col_header")
max(nchar(colnames(std_social_dosage_sumUnstd)))
#19
h5createDataset(file=fid,dataset="social_pruned_dosages/col_header/id",dims=dim(std_social_dosage_sumUnstd)[2],storage.mode='character',size=23)
h5write(obj=colnames(std_social_dosage_sumUnstd),file=fid,name="social_pruned_dosages/col_header/id")
h5createDataset(file=fid,dataset="social_pruned_dosages/matrix",dims=dim(std_social_dosage_sumUnstd),chunk = c(dim(std_social_dosage_sumUnstd)[1],ceiling(dim(std_social_dosage_sumUnstd)[2])/2))
h5write(obj=std_social_dosage_sumUnstd,file=fid,name="social_pruned_dosages/matrix")


load('data_dir/dosages/pruned_dosages/GRM.RData')
load('data_dir/simulations/power_analysis/fake_cage_cm.RData')
motch = match(rownames(fake_cage_cm),rownames(GRM))
any(is.na(motch))
# FALSE
GRM=GRM[motch,motch]

h5createGroup(fid,"GRM")
h5createGroup(fid,"GRM/pruned_dosages")

h5createGroup(fid,"GRM/pruned_dosages/row_header")
max(nchar(rownames(GRM)))
#15
h5createDataset(file=fid,dataset="GRM/pruned_dosages/row_header/sample_ID",dims=dim(GRM)[1],storage.mode='character',size=20)
h5write(obj=rownames(GRM),file=fid,name="GRM/pruned_dosages/row_header/sample_ID")

h5write(obj=GRM,file=fid,name="GRM/pruned_dosages/matrix")


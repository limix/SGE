my_load = function(arg) {
    load(paste('data_dir/simulations/simulations/simulations',arg,'.RData',sep=''))
    return(simulations)
}
res = lapply(1:67,FUN = my_load)

all_simulations = do.call('cbind',res)
all_simulations = all_simulations[,sample(1:dim(all_simulations)[2],replace = F)]
dim(all_simulations)
#[1]  1812 33500

library(rhdf5)
fid='data_dir/simulations/simulations/simulations.hdf5'
h5createFile(fid)

h5createGroup(fid,"simulations")

h5createDataset(file=fid,dataset="simulations/array",dims=dim(all_simulations))
h5write(obj=all_simulations,file=fid,name="simulations/array")

h5createGroup(fid,"simulations/rows_subjects")
max(nchar(rownames(all_simulations)))
#15
h5createDataset(file=fid,dataset="simulations/rows_subjects/outbred",dims=dim(all_simulations)[1],storage.mode='character',size=20)
h5write(obj=rownames(all_simulations),file=fid,name="simulations/rows_subjects/outbred")

#
h5createGroup(fid,"simulations/cols_measures")

max(nchar(colnames(all_simulations)))
#40
h5createDataset(file=fid,dataset="simulations/cols_measures/measures",dims=dim(all_simulations)[2],storage.mode='character',size=45)
h5write(obj=colnames(all_simulations),file=fid,name="simulations/cols_measures/measures")


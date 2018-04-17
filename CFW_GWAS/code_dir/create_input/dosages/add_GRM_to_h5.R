library(rhdf5)
fid='data_dir/CFWmice.h5'

load('data_dir/dosages/pruned_dosages/GRM.RData')

h5createGroup(fid,"GRM")
h5createGroup(fid,"GRM/pruned_dosages")

h5createGroup(fid,"GRM/pruned_dosages/row_header")
max(nchar(rownames(GRM)))
#15
h5createDataset(file=fid,dataset="GRM/pruned_dosages/row_header/sample_ID",dims=dim(GRM)[1],storage.mode='character',size=20)
h5write(obj=rownames(GRM),file=fid,name="GRM/pruned_dosages/row_header/sample_ID")

h5write(obj=GRM,file=fid,name="GRM/pruned_dosages/matrix")


library(rhdf5)
fid='data_dir/CFWmice.h5'
h5createFile(fid)

load('data_dir/phenotypes/all_pheno_data.RData')

######## phenotypes first
h5createGroup(fid,"phenotypes")
h5createGroup(fid,"phenotypes/row_header")
max(nchar(rownames(phenotypes_bc)))
#15
h5createDataset(file=fid,dataset="phenotypes/row_header/sample_ID",dims=dim(phenotypes_bc)[1],storage.mode='character',size=20)
h5write(obj=rownames(phenotypes_bc),file=fid,name="phenotypes/row_header/sample_ID")
max(nchar(cages),na.rm=T)
#7
h5createDataset(file=fid,dataset="phenotypes/row_header/cage",dims=length(cages),storage.mode='character',size=10)
h5write(obj=cages,file=fid,name="phenotypes/row_header/cage")
#h5createDataset(file=fid,dataset="phenotype/matrix",dims=dim(phenotypes))
h5write(obj=phenotypes_bc,file=fid,name="phenotypes/matrix")

h5createGroup(fid,"phenotypes/col_header")
max(nchar(colnames(phenotypes_bc)))
#41
h5createDataset(file=fid,dataset="phenotypes/col_header/phenotype_ID",dims=dim(phenotypes_bc)[2],storage.mode='character',size=45)
h5write(obj=colnames(phenotypes_bc),file=fid,name="phenotypes/col_header/phenotype_ID")

max(nchar(my_model_menu[,'my_covs']))
#56
h5createDataset(file=fid,dataset="phenotypes/col_header/covariatesUsed",dims=length(my_model_menu[,'my_covs']),storage.mode='character',size=60)
h5write(obj=my_model_menu[,'my_covs'],file=fid,name="phenotypes/col_header/covariatesUsed")

##### covariates now
h5createGroup(fid,"covariates")
h5createGroup(fid,"covariates/row_header")
max(nchar(rownames(cov_data)))
#15
h5createDataset(file=fid,dataset="covariates/row_header/sample_ID",dims=dim(cov_data)[1],storage.mode='character',size=20)
h5write(obj=rownames(cov_data),file=fid,name="covariates/row_header/sample_ID")

h5write(obj=cov_data,file=fid,name="covariates/matrix")

h5createGroup(fid,"covariates/col_header")
max(nchar(colnames(cov_data)))
#24
h5createDataset(file=fid,dataset="covariates/col_header/covariate_ID",dims=dim(cov_data)[2],storage.mode='character',size=30)
h5write(obj=colnames(cov_data),file=fid,name="covariates/col_header/covariate_ID")


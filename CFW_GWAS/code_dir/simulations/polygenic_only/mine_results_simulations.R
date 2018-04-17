library(parallel)

my_f=function(file) {
    VCs=try(read.delim(file,header=F,as.is=T))
    #if (grepl('no lines available in input',VCs)) file.remove(file)
    return(VCs)
}

dir = 'output_dir/VD_simulations/pruned_dosages_include_DGE_IGE_IEE_cageEffect/'
files=list.files(dir)
files=files[grep('.txt',files,fixed=T)]
length(files)

all_VCs=mclapply(paste(dir,'/',files,sep=''),my_f,mc.cores=20)
w=which(unlist(lapply(all_VCs,inherits,'try-error')))
length(w)
empty_files = files [w]
if (length(w)!=0) all_VCs=all_VCs[-w]

dims=do.call('rbind',lapply(all_VCs,FUN = dim))
w=which(dims[,1]!=50)
length(w)
dims[w,]

all_VCs=do.call('rbind',all_VCs)
dim(all_VCs)

colnames(all_VCs)=c('trait','sample_size','sample_size_cm','var_Ad','var_As','STE_Ad','STE_As','total_var','corr_Ads','STE_corr_Ads','corr_params','conv','LML','var_Ed','var_Es','corr_Eds','var_C') 
    
    
print(paste('All model converged: ',all(all_VCs[,'conv']=='True'),sep=''))
w=which(all_VCs[,'conv']=='False')
length(w)
#0
if (length(w)!=0) all_VCs=all_VCs[-w,]
    
for (col in c('var_Ad','var_As','corr_Ads','var_Ed','var_Es','corr_Eds','var_C')) {
    all_VCs[which(all_VCs[,col]==(-999)),col]=0
}

all_VCs[,c('var_Ad','var_As','var_Ed','var_Es','var_C')] = all_VCs[,c('var_Ad','var_As','var_Ed','var_Es','var_C')] / all_VCs[,'total_var'] * 100

save(all_VCs, file = 'output_dir/VD_simulations/pruned_dosages_include_DGE_IGE_IEE_cageEffect/all_VCs.RData')

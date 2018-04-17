#path='output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/w_STE'
#constrained = NULL

#to calculate significance of aggregate SGE
#path = 'output_dir/VD/pruned_dosages_include_DGE_IEE_cageEffect'
#constrained = NULL

# significance of aggregate DGE
path = 'output_dir/VD/pruned_dosages_include_IGE_IEE_cageEffect'
constrained = NULL


#to calculate significance of rho != 0
#path = 'output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/constrained1'
#constrained = 1

#for QQ plots
#path='output_dir/VD/pruned_dosages_include_DGE_cageEffect'
#constrained = NULL

#for QQ plots
#path='output_dir/VD/pruned_dosages_include_cageEffect'
#constrained = NULL


files=list.files(path)
if (!is.null(constrained)) files=files[grep(paste('_constrained',constrained,'.txt',sep=''),files,fixed=T)] else files=files[grepl('.txt',files,fixed=T) & !grepl('constrained',files,fixed=T)]
if (length(files)==0) {print('no files') ; next}
all_VCs=NULL
for (file in files) {
    VCs=try(read.delim(paste(path,'/',file,sep=''),header=F,as.is=T))
    if (inherits(VCs,'try-error')) next #stop('1') 
    all_VCs=rbind(all_VCs,VCs)
}

    if (is.null(constrained)) colnames(all_VCs)=c('trait','sample_size','sample_size_cm','var_Ad','var_As','STE_Ad','STE_As','total_var','corr_Ads','STE_corr_Ads','corr_params','conv','LML','var_Ed','var_Es','corr_Eds','var_C')  else colnames(all_VCs)=c('trait','sample_size','sample_size_cm','var_Ad','var_As','total_var','corr_Ads','conv','LML')

if (is.null(constrained)) rownames(all_VCs)=all_VCs[,'trait']
for (col in colnames(all_VCs)) {
    all_VCs[which(all_VCs[,col]==(-999)),col]=NA
}
    
print(paste('All model converged: ',all(all_VCs[,'conv']=='True'),sep=''))
    
if (is.null(constrained)) all_VCs[,c('var_Ad','var_As','var_Ed','var_Es','var_C')] = all_VCs[,c('var_Ad','var_As','var_Ed','var_Es','var_C')] / all_VCs[,'total_var'] else all_VCs[,c('var_Ad','var_As')] = all_VCs[,c('var_Ad','var_As')] / all_VCs[,'total_var']


phenos_to_run = read.delim('data_dir/phenotypes/colnames_phenotypes_bc_to_run.txt',as.is =T, header = F)[,1]
motch = match(phenos_to_run, all_VCs[,'trait'])
if (any(is.na(motch))) stop('pb')
w = which(all_VCs[,'trait'] %in% phenos_to_run)
all_VCs = all_VCs[w,]

if (is.null(constrained)) {
    all_VCs=all_VCs[order(all_VCs[,'var_As'],decreasing=T),]
    save(all_VCs,file=paste(path,'/all_VCs.RData',sep=''))
} else {
    save(all_VCs,file=paste(path,'/all_VCs_constrained',constrained,'.RData',sep=''))
}
print(dim(all_VCs)[1])




# using output_dir/pruned_dosages_include_DGE_IGE_IEE_cageEffect/all_VCs.RData, compare DGE and SGE of phenotypes with and without body weight correction
g = grep('_BWcorr',rownames(all_VCs))
BWcorr = rownames(all_VCs)[g]
no_BWcorr = sub('_BWcorr','',BWcorr)
    
pdf('plots_dir/check_collider.pdf')    
par(mfrow = c(1,2))
plot(all_VCs[no_BWcorr,'var_Ad'],all_VCs[BWcorr,'var_Ad'])
abline(0,1,col = 'red')
plot(all_VCs[no_BWcorr,'var_As'],all_VCs[BWcorr,'var_As'])
abline(0,1,col = 'red')
dev.off()

# no collider effect at genome-wide level
# so will look at BW_corr measures only




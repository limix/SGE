load('data_dir/phenotypes/my_model_menu.RData')

my_model_menu = my_model_menu[order(my_model_menu[,'Name']),]

imputed_cage_info = read.delim('data_dir/Supplementary_Table3.xls',as.is = T, header = T)
motch=match(rownames(phenotypes),rownames(imputed_cage_info))
any(is.na(motch))
#FALSE
imputed_cage_info = imputed_cage_info[motch,]
cages = imputed_cage_info[,1]
names(cages) = rownames(imputed_cage_info)


############ boxcox
data = cbind(phenotypes,covariates[,!colnames(covariates) %in% colnames(phenotypes)])
library(lme4)
library(parallel)
library(MASS)
lambda_limit=2
for (i in 1:dim(my_model_menu)[1]) {
	pheno_name=my_model_menu[i,'Name']
	if (!pheno_name %in% colnames(data)) {print(paste(pheno_name,'not in colnames(data)')) ; stop()}
	m=min(data[,pheno_name],na.rm=T)
	if (m<=0) data[,pheno_name]=data[,pheno_name]+abs(m)+1
	
	these_covs = strsplit(my_model_menu[pheno_name,'my_covs'],',')[[1]]
	form=as.formula(paste(pheno_name,paste('~',paste(these_covs, collapse = '+'))))
	#print(form)
	bc=boxcox(form,data=data,plotit=F,lambda=seq(-lambda_limit,lambda_limit,by=0.1))
	lambda=bc$x[which.max(bc$y)]
	#retain BMC.osteoporosis as it is a binary phenotype (!! only very large QTLs will be considered as few rats with osteoporosis)
	if (abs(lambda) == lambda_limit & pheno_name!='BMC.osteoporosis') {print(paste(pheno_name,'excluded because cannot be transformed to gaussian')) ; data[,pheno_name] = NA} else {
		if (lambda==0) data[,pheno_name]=log(data[,pheno_name]) else 
			data[,pheno_name]=(data[,pheno_name]^lambda-1)/lambda
		if (any(abs(data[,pheno_name]-(-999))<3,na.rm=T)) {print('pb') ; stop('pb')}
	}
	data[is.na(data[,pheno_name]),pheno_name]=(-999)
}
#[1] "Bioch.Sodium excluded because cannot be transformed to gaussian"
#[1] "Bioch.Sodium_BWcorr excluded because cannot be transformed to gaussian"
#[1] "BMC.Max excluded because cannot be transformed to gaussian"
#[1] "BMC.Mode excluded because cannot be transformed to gaussian"
#[1] "BMC.Skew excluded because cannot be transformed to gaussian"
#[1] "BMC.Width excluded because cannot be transformed to gaussian"
#[1] "Cardio.ECG.Heart_Rate excluded because cannot be transformed to gaussian"
#[1] "Cardio.ECG.Heart_Rate_BWcorr excluded because cannot be transformed to gaussian"
#[1] "Diss.Body.Length excluded because cannot be transformed to gaussian"
#[1] "Diss.Body.Length_BWcorr excluded because cannot be transformed to gaussian"
#[1] "Haem.abs_basos excluded because cannot be transformed to gaussian"
#[1] "Haem.abs_eos excluded because cannot be transformed to gaussian"
#[1] "Haem.abs_lucs excluded because cannot be transformed to gaussian"
#[1] "Haem.abs_mono excluded because cannot be transformed to gaussian"
#[1] "Haem.BASO_percent excluded because cannot be transformed to gaussian"
#[1] "Haem.HCT excluded because cannot be transformed to gaussian"
#[1] "Haem.LYM_percent excluded because cannot be transformed to gaussian"
#[1] "Haem.measHGB excluded because cannot be transformed to gaussian"
#[1] "Haem.RBC excluded because cannot be transformed to gaussian"
#[1] "OFT.Periphery.Time excluded because cannot be transformed to gaussian"
#[1] "Sleep.long_sleep excluded because cannot be transformed to gaussian"
#[1] "Sleep.long_sleep_BWcorr excluded because cannot be transformed to gaussian"
#[1] "Sleep.Percent_wake_over_17min excluded because cannot be transformed to gaussian"
#[1] "Sleep.Percent_wake_over_17min_BWcorr excluded because cannot be transformed to gaussian"
#[1] "SPPI.ln_pc excluded because cannot be transformed to gaussian"
#[1] "SPPI.ln_pc_BWcorr excluded because cannot be transformed to gaussian"

my_all = function(col) {
	return(all(col == (-999)))
}
w = which(apply(data,FUN = my_all, MAR = 2))
sort(colnames(data)[w])
# corresponds to what was printed above
data = data[,-w]
# intersect because data also have covariates and phenotypes has measures that failed boxcox
keep = intersect(colnames(phenotypes), colnames(data))
#data has transformed data
phenotypes_bc = data[,keep]

all (rownames(phenotypes_bc) == rownames(phenotypes))
#TRUE

motch = match(colnames(phenotypes_bc),my_model_menu[,'Name'])
any(is.na(motch))
#[1] FALSE
my_model_menu = my_model_menu[motch,]
all (colnames(phenotypes_bc) == my_model_menu[,'Name'])
#[1] TRUE

pdf('plots_dir/phenotypes_boxcoxedVSraw.pdf',width=15)
par(mfrow=c(1,3))
for (measure in colnames(phenotypes_bc)) {
	w = which(phenotypes_bc[,measure]!=(-999))
	truehist(phenotypes[w,measure])
	truehist(phenotypes_bc[w,measure])
	plot(phenotypes[w,measure],phenotypes_bc[w,measure],main = measure,sub = paste(sum(is.na(phenotypes[,measure])),sum(phenotypes_bc[,measure] == (-999))))
}
dev.off()


############ get numerical covariates matrix
my_strsplit = function(string) {
	return(strsplit(string,',')[[1]])
}
unique_covs = unique(unlist(sapply(my_model_menu[,'my_covs'],FUN = my_strsplit)))
unique_covs
# [1] "Sex"                      "cage_density"            
# [3] "Weight.Hypo"              "Weight.Diss"             
# [5] "BMC.osteoporosis"         "Weight.Startle"          
# [7] "FC.Box.Cue"               "Hypoxia.Chamber"         
# [9] "Hypoxia.Time.Hour.Minute" "Tibia.Length"            
#[11] "Neuro.Nb_slices_DCX"      "Neuro.Nb_slices_Ki67"    

covariates[,'Sex'] = as.numeric(as.factor(covariates[,'Sex']))
#M recoded as 2 and females as 1
cov_data = covariates[,c('Sex',"Weight.Hypo","Weight.Diss",'BMC.osteoporosis',"Weight.Startle" ,'Hypoxia.Time.Hour.Minute','Tibia.Length','Neuro.Nb_slices_DCX','Neuro.Nb_slices_Ki67')]
for (covariate in c('Hypoxia.Chamber','FC.Box.Cue')) {
	for (val in na.exclude(unique(covariates[,covariate]))[-1]) {
    	cov=matrix(0,,nrow = dim(covariates)[1],ncol=1)
    	cov[which(covariates[,covariate]==val),1]=1
    	colnames(cov) = covariate
    	cov_data=cbind(cov_data,cov)
	}
}

for(i in 1:dim(cov_data)[2]) {
	if (any(abs(cov_data[,i]-(-999))<3,na.rm=T)) stop('pb')
    cov_data[is.na(cov_data[,i]),i]=(-999)
}
cov_data=as.matrix(cov_data)

phenotypes_bc = as.matrix(phenotypes_bc)
str(phenotypes_bc)
#num
str(cov_data)
#num

save(cages,phenotypes_bc,cov_data,my_model_menu,file = 'data_dir/phenotypes/all_pheno_data.RData')




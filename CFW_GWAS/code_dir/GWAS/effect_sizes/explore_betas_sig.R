library(rhdf5)

load('output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/w_STE/all_VCs.RData')

#leave as is: first part is for variance explained, second part (below) for betas
type = 'var_expl'

DGE_dir = paste('output_dir/effect_sizes_additive_model/betas_conditioning/',type,'/DGE/',sep='')
DGE_betas_files = list.files(DGE_dir)
DGE_betas_files = DGE_betas_files[grep('txt',DGE_betas_files)]
SGE_dir = paste('output_dir/effect_sizes_additive_model/betas_conditioning/',type,'/SGE/',sep='')
SGE_betas_files = list.files(SGE_dir)
SGE_betas_files = SGE_betas_files[grep('txt',SGE_betas_files)]
inter_betas = intersect(DGE_betas_files,SGE_betas_files)
length(inter_betas)
#72

#all_res will have total var explained by all QTLs for a phenotype
all_res_DGE = NULL
all_res_SGE = NULL
#indiv will have var explained by individual QTLs
indiv_DGE_QTLs_betas = NULL
indiv_SGE_QTLs_betas = NULL
nomes_DGE = c()
nomes_SGE = c()
for (pheno_name in sub('.txt','',inter_betas)) {

	total_var_null = all_VCs[pheno_name,'total_var']
	DGE_herit = all_VCs[pheno_name,'var_Ad']
	SGE_herit = all_VCs[pheno_name,'var_As']

	res_DGE_QTLs = read.table(paste(DGE_dir,pheno_name,'.txt',sep=''),as.is = T)
	#position id top QTL
	w = grep('direct',res_DGE_QTLs)
	if (length(w) >= 1) {
		nb_non_QTL_covs = (dim(res_DGE_QTLs)[2] -3 -2*length(w) -5)/2
		nb_non_QTL_covs_before = w[1] - 3 - 1
		nb_non_QTL_covs_after = nb_non_QTL_covs - nb_non_QTL_covs_before
		w_mean = grep('mean',res_DGE_QTLs)

		if (nb_non_QTL_covs_after!=0) colnames(res_DGE_QTLs) = c('trait','sample','sample_cm',paste('cov',1:nb_non_QTL_covs_before,'_id',sep=''),paste('sig_id',1:length(w),sep=''),paste('cov',(nb_non_QTL_covs_before+1):(nb_non_QTL_covs_before+nb_non_QTL_covs_after),'_id',sep=''),paste('cov',1:nb_non_QTL_covs_before,'_beta',sep=''),paste('sig_beta',1:length(w),sep=''),paste('cov',(nb_non_QTL_covs_before+1):(nb_non_QTL_covs_before+nb_non_QTL_covs_after),'_beta',sep=''),'DGE_herit','SGE_herit','total_var','conv','LML') else {
			colnames(res_DGE_QTLs) = c('trait','sample','sample_cm',paste('cov',1:nb_non_QTL_covs_before,'_id',sep=''),paste('sig_id',1:length(w),sep=''),paste('cov',1:nb_non_QTL_covs_before,'_beta',sep=''),paste('sig_beta',1:length(w),sep=''),'DGE_herit','SGE_herit','total_var','conv','LML')
		}
		colnames(res_DGE_QTLs)[w_mean] = 'mean_id'
		colnames(res_DGE_QTLs)[w_mean + nb_non_QTL_covs + length(w)] = 'mean_intercept'

		betas = unlist(res_DGE_QTLs[1,grep('beta',colnames(res_DGE_QTLs))])
		#genotypes standardized to var 1 when type = 'var_expl' so variance explained by QTL is beta^2
		this_total_var = unlist(res_DGE_QTLs[1,'total_var']) + sum(betas^2)
		#proportion of variance explained
		betas = betas^2 / this_total_var
		beta_sig_DGE = betas[grep('sig_beta',names(betas))]
		for (i in 1:length(beta_sig_DGE)) {
			indiv_DGE_QTLs_betas = rbind(indiv_DGE_QTLs_betas,c(pheno_name,res_DGE_QTLs[,paste('sig_id',i,sep='')],beta_sig_DGE[paste('sig_beta',i,sep='')]))
		}
		beta_sig_DGE = sum(beta_sig_DGE)
		all_res_DGE = rbind(all_res_DGE,c(beta_sig_DGE,DGE_herit,this_total_var))
		nomes_DGE = c(nomes_DGE,pheno_name)
	}

	res_SGE_QTLs = read.table(paste(SGE_dir,pheno_name,'.txt',sep=''),as.is = T)
	#position id top QTL
	w = grep('social',res_SGE_QTLs)
	if (length(w) >= 1) {
		nb_non_QTL_covs = (dim(res_SGE_QTLs)[2] -3 -2*length(w) -5)/2
		nb_non_QTL_covs_before = w[1] - 3 - 1
		nb_non_QTL_covs_after = nb_non_QTL_covs - nb_non_QTL_covs_before
		w_mean = grep('mean',res_SGE_QTLs)

		if (nb_non_QTL_covs_after!=0) colnames(res_SGE_QTLs) = c('trait','sample','sample_cm',paste('cov',1:nb_non_QTL_covs_before,'_id',sep=''),paste('sig_id',1:length(w),sep=''),paste('cov',(nb_non_QTL_covs_before+1):(nb_non_QTL_covs_before+nb_non_QTL_covs_after),'_id',sep=''),paste('cov',1:nb_non_QTL_covs_before,'_beta',sep=''),paste('sig_beta',1:length(w),sep=''),paste('cov',(nb_non_QTL_covs_before+1):(nb_non_QTL_covs_before+nb_non_QTL_covs_after),'_beta',sep=''),'DGE_herit','SGE_herit','total_var','conv','LML') else {
			colnames(res_SGE_QTLs) = c('trait','sample','sample_cm',paste('cov',1:nb_non_QTL_covs_before,'_id',sep=''),paste('sig_id',1:length(w),sep=''),paste('cov',1:nb_non_QTL_covs_before,'_beta',sep=''),paste('sig_beta',1:length(w),sep=''),'DGE_herit','SGE_herit','total_var','conv','LML')
		}
		colnames(res_SGE_QTLs)[w_mean] = 'mean_id'
		colnames(res_SGE_QTLs)[w_mean + nb_non_QTL_covs + length(w)] = 'mean_intercept'

		betas = unlist(res_SGE_QTLs[1,grep('beta',colnames(res_SGE_QTLs))])
		#genotypes standardized to var 1 when type = 'var_expl' so variance explained by QTL is beta^2
		this_total_var = unlist(res_SGE_QTLs[1,'total_var']) + sum(betas^2)
		#proportion of variance explained
		betas = betas^2 / this_total_var
		beta_sig_SGE = betas[grep('sig_beta',names(betas))]
		for (i in 1:length(beta_sig_SGE)) {
			indiv_SGE_QTLs_betas = rbind(indiv_SGE_QTLs_betas,c(pheno_name,res_SGE_QTLs[,paste('sig_id',i,sep='')],beta_sig_SGE[paste('sig_beta',i,sep='')]))
		}
		beta_sig_SGE = sum(beta_sig_SGE)
		all_res_SGE = rbind(all_res_SGE,c(beta_sig_SGE,SGE_herit,this_total_var))
		nomes_SGE = c(nomes_SGE,pheno_name)
	}	
}
colnames(all_res_DGE) = c('var_expl_DGE_QTLs','herit_DGE','total_var')
rownames(all_res_DGE) = nomes_DGE
all_res_DGE = as.data.frame(all_res_DGE,stringsAsFactors =F)
all_res_DGE[,'var_expl_DGE_QTLs'] = as.numeric(all_res_DGE[,'var_expl_DGE_QTLs'])
all_res_DGE[,'herit_DGE'] = as.numeric(all_res_DGE[,'herit_DGE'])
all_res_DGE[,'total_var'] = as.numeric(all_res_DGE[,'total_var'])

colnames(all_res_SGE) = c('var_expl_SGE_QTLs','herit_SGE','total_var')
rownames(all_res_SGE) = nomes_SGE
all_res_SGE = as.data.frame(all_res_SGE,stringsAsFactors =F)
all_res_SGE[,'var_expl_SGE_QTLs'] = as.numeric(all_res_SGE[,'var_expl_SGE_QTLs'])
all_res_SGE[,'herit_SGE'] = as.numeric(all_res_SGE[,'herit_SGE'])
all_res_SGE[,'total_var'] = as.numeric(all_res_SGE[,'total_var'])

indiv_DGE_QTLs_betas = as.data.frame(indiv_DGE_QTLs_betas, stringsAsFactors = F)
indiv_DGE_QTLs_betas[,'sig_beta1'] = as.numeric(indiv_DGE_QTLs_betas[,'sig_beta1'])
dim(indiv_DGE_QTLs_betas)
#[1] 118   3
mean(indiv_DGE_QTLs_betas[,'sig_beta1'])
#0.02750621 var expl

indiv_SGE_QTLs_betas = as.data.frame(indiv_SGE_QTLs_betas, stringsAsFactors = F)
indiv_SGE_QTLs_betas[,'sig_beta1'] = as.numeric(indiv_SGE_QTLs_betas[,'sig_beta1'])
dim(indiv_SGE_QTLs_betas)
#[1] 21  3
mean(indiv_SGE_QTLs_betas[,'sig_beta1'])
#0.01917476 var expl

indiv_DGE_QTLs_betas[,'sig_beta1'] = indiv_DGE_QTLs_betas[,'sig_beta1'] *100
indiv_SGE_QTLs_betas[,'sig_beta1'] = indiv_SGE_QTLs_betas[,'sig_beta1'] *100

rx = range(c(indiv_DGE_QTLs_betas[,'sig_beta1'],indiv_SGE_QTLs_betas[,'sig_beta1']))
brokes = seq(0,rx[2]+5,by = 5)

save(indiv_DGE_QTLs_betas,indiv_SGE_QTLs_betas,all_res_DGE,all_res_SGE,file = paste('output_dir/effect_sizes_additive_model/betas_conditioning/',type,'/betas.RData',sep=''))


####

# will actually look at betas squared
type = 'allelic_effect'
DGE_dir = paste('output_dir/effect_sizes_additive_model/betas_conditioning/',type,'/DGE/',sep='')
SGE_dir = paste('output_dir/effect_sizes_additive_model/betas_conditioning/',type,'/SGE/',sep='')
#load(paste('output_dir/effect_sizes/betas/var_expl/betas.RData',sep=''))
load(paste('output_dir/effect_sizes_additive_model/betas_conditioning/var_expl/betas.RData',sep=''))
indiv_DGE_QTLs_betas = NULL
indiv_SGE_QTLs_betas = NULL
for (pheno_name in rownames(all_res_DGE)) {

	res_DGE_QTLs = read.table(paste(DGE_dir,pheno_name,'.txt',sep=''),as.is = T)
	#position id top QTL
	w = grep('direct',res_DGE_QTLs)
	if (length(w) >= 1) {
		nb_non_QTL_covs = (dim(res_DGE_QTLs)[2] -3 -2*length(w) -5)/2
		nb_non_QTL_covs_before = w[1] - 3 - 1
		nb_non_QTL_covs_after = nb_non_QTL_covs - nb_non_QTL_covs_before
		w_mean = grep('mean',res_DGE_QTLs)

		if (nb_non_QTL_covs_after!=0) colnames(res_DGE_QTLs) = c('trait','sample','sample_cm',paste('cov',1:nb_non_QTL_covs_before,'_id',sep=''),paste('sig_id',1:length(w),sep=''),paste('cov',(nb_non_QTL_covs_before+1):(nb_non_QTL_covs_before+nb_non_QTL_covs_after),'_id',sep=''),paste('cov',1:nb_non_QTL_covs_before,'_beta',sep=''),paste('sig_beta',1:length(w),sep=''),paste('cov',(nb_non_QTL_covs_before+1):(nb_non_QTL_covs_before+nb_non_QTL_covs_after),'_beta',sep=''),'DGE_herit','SGE_herit','total_var','conv','LML') else {
			colnames(res_DGE_QTLs) = c('trait','sample','sample_cm',paste('cov',1:nb_non_QTL_covs_before,'_id',sep=''),paste('sig_id',1:length(w),sep=''),paste('cov',1:nb_non_QTL_covs_before,'_beta',sep=''),paste('sig_beta',1:length(w),sep=''),'DGE_herit','SGE_herit','total_var','conv','LML')
		}
		colnames(res_DGE_QTLs)[w_mean] = 'mean_id'
		colnames(res_DGE_QTLs)[w_mean + nb_non_QTL_covs + length(w)] = 'mean_intercept'

		betas = unlist(res_DGE_QTLs[1,grep('beta',colnames(res_DGE_QTLs))])
		#absolute additive effect of the ref alleles on the phenotype in standard deviations (https://www.nature.com/articles/nrg.2017.101#f4)
		
		betas = betas^2
		# now something more readily usable for comparison of power between sge and dgeGWAS		
		this_total_var = all_res_DGE[pheno_name,'total_var']
		betas = abs(betas) / this_total_var
		beta_sig_DGE = betas[grep('sig_beta',names(betas))]
		for (i in 1:length(beta_sig_DGE)) {
			indiv_DGE_QTLs_betas = rbind(indiv_DGE_QTLs_betas,c(pheno_name,res_DGE_QTLs[,paste('sig_id',i,sep='')],beta_sig_DGE[paste('sig_beta',i,sep='')]))
		}
	}
}
for (pheno_name in rownames(all_res_SGE)) {
	res_SGE_QTLs = read.table(paste(SGE_dir,pheno_name,'.txt',sep=''),as.is = T)
	#position id top QTL
	w = grep('social',res_SGE_QTLs)
	if (length(w) >= 1) {
		nb_non_QTL_covs = (dim(res_SGE_QTLs)[2] -3 -2*length(w) -5)/2
		nb_non_QTL_covs_before = w[1] - 3 - 1
		nb_non_QTL_covs_after = nb_non_QTL_covs - nb_non_QTL_covs_before
		w_mean = grep('mean',res_SGE_QTLs)

		if (nb_non_QTL_covs_after!=0) colnames(res_SGE_QTLs) = c('trait','sample','sample_cm',paste('cov',1:nb_non_QTL_covs_before,'_id',sep=''),paste('sig_id',1:length(w),sep=''),paste('cov',(nb_non_QTL_covs_before+1):(nb_non_QTL_covs_before+nb_non_QTL_covs_after),'_id',sep=''),paste('cov',1:nb_non_QTL_covs_before,'_beta',sep=''),paste('sig_beta',1:length(w),sep=''),paste('cov',(nb_non_QTL_covs_before+1):(nb_non_QTL_covs_before+nb_non_QTL_covs_after),'_beta',sep=''),'DGE_herit','SGE_herit','total_var','conv','LML') else {
			colnames(res_SGE_QTLs) = c('trait','sample','sample_cm',paste('cov',1:nb_non_QTL_covs_before,'_id',sep=''),paste('sig_id',1:length(w),sep=''),paste('cov',1:nb_non_QTL_covs_before,'_beta',sep=''),paste('sig_beta',1:length(w),sep=''),'DGE_herit','SGE_herit','total_var','conv','LML')
		}
		colnames(res_SGE_QTLs)[w_mean] = 'mean_id'
		colnames(res_SGE_QTLs)[w_mean + nb_non_QTL_covs + length(w)] = 'mean_intercept'

		betas = unlist(res_SGE_QTLs[1,grep('beta',colnames(res_SGE_QTLs))])
		betas = betas^2
		this_total_var = all_res_SGE[pheno_name,'total_var']
		betas = abs(betas) / this_total_var
		beta_sig_SGE = betas[grep('sig_beta',names(betas))]
		for (i in 1:length(beta_sig_SGE)) {
			indiv_SGE_QTLs_betas = rbind(indiv_SGE_QTLs_betas,c(pheno_name,res_SGE_QTLs[,paste('sig_id',i,sep='')],beta_sig_SGE[paste('sig_beta',i,sep='')]))
		}
	}	
}
indiv_DGE_QTLs_betas = as.data.frame(indiv_DGE_QTLs_betas, stringsAsFactors = F)
indiv_DGE_QTLs_betas[,'sig_beta1'] = as.numeric(indiv_DGE_QTLs_betas[,'sig_beta1'])
indiv_SGE_QTLs_betas = as.data.frame(indiv_SGE_QTLs_betas, stringsAsFactors = F)
indiv_SGE_QTLs_betas[,'sig_beta1'] = as.numeric(indiv_SGE_QTLs_betas[,'sig_beta1'])

save(indiv_DGE_QTLs_betas,indiv_SGE_QTLs_betas,file = paste('output_dir/effect_sizes_additive_model/betas_conditioning/',type,'/betas.RData',sep=''))

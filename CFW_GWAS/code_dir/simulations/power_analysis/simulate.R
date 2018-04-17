#SNPsCanUse13_colMax.RData to SNPsCanUse24_colMax.RData
#genotypic_correlations_gs3_colMax.RData

options(warn = 2)
load('data_dir/simulations/power_analysis/MAFs.RData')


# structure of sims_design below (ie which parameters will be used in simulations) needs to match hard-coded parts of efficient_map_SGE_single_snp.py and draw_sims_wPolygenic.R

sims_design = matrix(nrow = 0, ncol = 5)
colnames(sims_design) = c('gs','MAF','corsP','sgeQTL_a','dgeQTL_a')
row_sims_design = 0
for (gs in 2:6) {
	print(gs)
	if (gs>2) rm(list = c('low_corsP','high_corsP'))
	load(paste('data_dir/simulations/power_analysis/genotypic_correlations_gs',gs,'.RData', sep=''))
	for (maf in c('low','mid','high')) {
		for (corsP in c('low','high')) {
			SNPsCanUse = intersect(get(paste(maf,'MAFs',sep='_')),get(paste(corsP,'corsP',sep='_')))
			print(length(SNPsCanUse))

			#no DGE QTL
			sims_design = rbind(sims_design,c(gs,maf,corsP,0.2,0))
			row_sims_design = row_sims_design + 1
			save(SNPsCanUse,file = paste('data_dir/simulations/power_analysis/SNPsCanUse/SNPsCanUse',row_sims_design,'.RData', sep=''))

			#lareg DGE QTL
			sims_design = rbind(sims_design,c(gs,maf,corsP,0.2,1))
			row_sims_design = row_sims_design + 1
			#yes same SNPs saved twice
			save(SNPsCanUse,file = paste('data_dir/simulations/power_analysis/SNPsCanUse/SNPsCanUse',row_sims_design,'.RData', sep=''))
		}
	}
}
sims_design = as.data.frame(sims_design,stringsAsFactors = F)
sims_design[,'gs'] = as.numeric(sims_design[,'gs'])
sims_design[,'sgeQTL_a'] = as.numeric(sims_design[,'sgeQTL_a'])
sims_design[,'dgeQTL_a'] = as.numeric(sims_design[,'dgeQTL_a'])
save(sims_design, file = 'data_dir/simulations/power_analysis/sims_design.RData')

dim(sims_design)

#note: I have added in other iterations 12 rows for additive SGE model (gs 3) and 12 rows for colMax SGE model (not in paper)

#always at least 1,000 SNPsCanUse


### now if want to get input for simulations with polygenic effects, will need different input for diff group size (as cages different)

#run for gs in 2,3,4,5,6
gs = 3

library(MASS)
library(rhdf5)

load('data_dir/dosages/pruned_dosages/GRM.RData')
load('data_dir/simulations/power_analysis/fake_cage_cm.RData')

motch = match(rownames(fake_cage_cm),rownames(GRM))
any(is.na(motch))
# FALSE
GRM=GRM[motch,motch]

#genetic matrices first
#K
n=dim(GRM)[1]
vecti=matrix(1,ncol=1,nrow=n)
p=diag(n)-(vecti%*%t(vecti))/n
num=(n-1)
denom=sum(diag(p%*%GRM%*%p))
my_scaling_factor_K=num/denom
scaled_K = my_scaling_factor_K * GRM

#ZKZ'
Z=matrix(0,nrow=length(fake_cage_cm[,paste('cage_gs',gs,sep='')]),ncol=length(fake_cage_cm[,paste('cage_gs',gs,sep='')]))
for (i in 1:length(fake_cage_cm[,paste('cage_gs',gs,sep='')])) {
    that_cage=fake_cage_cm[,paste('cage_gs',gs,sep='')][i]
    w=which(fake_cage_cm[,paste('cage_gs',gs,sep='')]==that_cage)
    w=w[w!=i]
    if (length(w)!=(gs - 1)) stop('pb')
    Z[i,w]=1
}
ZKZt= (Z %*% GRM) %*% t(Z)
n=dim(ZKZt)[1]
vecti=matrix(1,ncol=1,nrow=n)
p=diag(n)-(vecti%*%t(vecti))/n
num=(n-1)
denom=sum(diag(p%*%ZKZt%*%p))
my_scaling_factor_ZKZt=num/denom
scaled_K_cm = my_scaling_factor_ZKZt * GRM

#cross cov
my_scaling_factor_crossK = sqrt(my_scaling_factor_K*my_scaling_factor_ZKZt)
scaled_crossK = my_scaling_factor_crossK * GRM

#environmental matrices next
#Ed already scaled
env = diag(length(fake_cage_cm[,paste('cage_gs',gs,sep='')]))
my_scaling_factor_env=1
scaled_env = my_scaling_factor_env * env

#Es
ZZt= (Z %*% t(Z))
n=dim(ZZt)[1]
vecti=matrix(1,ncol=1,nrow=n)
p=diag(n)-(vecti%*%t(vecti))/n
num=(n-1)
denom=sum(diag(p%*%ZZt%*%p))
my_scaling_factor_ZZt=num/denom
scaled_env_cm=my_scaling_factor_ZZt * env

#cross env cov
my_scaling_factor_crossE = sqrt(my_scaling_factor_env*my_scaling_factor_ZZt)
scaled_env_cross = my_scaling_factor_crossE*env

#cage matrix finally (corresponds to WW in python code I think)
cage_mat = matrix(0,dim(GRM)[1],dim(GRM)[1])
for (cage in unique(fake_cage_cm[,paste('cage_gs',gs,sep='')])) {
    w=which(fake_cage_cm[,paste('cage_gs',gs,sep='')]==cage)
    cage_mat[w,w]=1
}
all((diag(ZZt)+1) == apply(cage_mat,FUN=sum,MAR=1))
#TRUE
#both ZZt and cage_mat are 1911 x 1911
#which(cage_mat[1,]==1)
#[1]   1  23 110 176 277
#which(cage_mat[23,]==1)
#[1]   1  23 110 176 277
#ZZt[1,1]
#[1] 4
#ZZt[1,23]
#[1] 3
#ZZt[1,110]
#[1] 3
#ZZt[1,176]
#[1] 3
#ZZt[1,277]
#[1] 3
#rest is 0
#ZZt is # of cagemates

n=dim(cage_mat)[1]
vecti=matrix(1,ncol=1,nrow=n)
p=diag(n)-(vecti%*%t(vecti))/n
num=(n-1)
denom=sum(diag(p%*%cage_mat%*%p))
my_scaling_factor_cage_mat=num/denom
scaled_cage_mat=my_scaling_factor_cage_mat * cage_mat

save(GRM, Z, scaled_K, scaled_K_cm, scaled_crossK, scaled_env, scaled_env_cm, scaled_env_cross, scaled_cage_mat, file = paste('data_dir/simulations/power_analysis/image_simulations_gs',gs,'.RData',sep=''))


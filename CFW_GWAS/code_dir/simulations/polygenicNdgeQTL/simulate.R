library(MASS)
library(rhdf5)

load('data_dir/dosages/pruned_dosages/GRM.RData')
load('data_dir/phenotypes/imputed_cage_info.RData')
load('data_dir/dosages/close_pairs_GRM.RData')

inter=intersect(colnames(GRM),rownames(imputed_cage_info))
inter = inter[!inter %in% all_cagemates_close_pairs]
length(inter)
#[1] 1812
GRM=GRM[match(inter,rownames(GRM)),match(inter,colnames(GRM))]
imputed_cage_info=imputed_cage_info[match(inter,rownames(imputed_cage_info)),]
cage_cm = imputed_cage_info[,'cage_WH']

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
Z=matrix(0,nrow=length(cage_cm),ncol=length(cage_cm))
for (i in 1:length(cage_cm)) {
    that_cage=cage_cm[i]
    w=which(cage_cm==that_cage)
    w=w[w!=i]
    if (length(w)!=2) stop('pb')
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
env = diag(length(cage_cm))
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
for (cage in unique(cage_cm)) {
    w=which(cage_cm==cage)
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

load('data_dir/dosages/pruned_dosages/unstd_my_final_dosages.RData')

#for info only
load("output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/all_VCs.RData")
medians = apply(all_VCs[,-c(1,12)],FUN=median,MAR=2)
round(medians['var_Ad']*100,digit=0)
#15 vs 9 for all
round(medians['var_As']*100,digit=0)
#8 vs 1
round(medians['corr_Ads'],digit=2)
#0.47 vs 0.44
round(medians['var_Ed'],digit=2)
#0.22 vs 0.3
round(medians['var_Es'],digit=2)
#0.16 vs 0.21
round(medians['corr_Eds'],digit=2)
#-0.97 vs -0.95
round(medians['var_C'],digit=2)
#26 vs 32

#var_Ad, var_As, corr_Ads, var_Ed, var_Es, corr_Eds, var_Cage
high_polygenic = c(20,20,0.5,15,15,-0.97,25)
low_polygenic = c(5,5,0.5,30,30,-0.97,25)
vars_dgeQTL = c(0,5,20,50)

liste_effects = NULL
for (var_dgeQTL in vars_dgeQTL) {
	liste_effects = rbind(liste_effects,c(low_polygenic,var_dgeQTL))
	liste_effects = rbind(liste_effects,c(high_polygenic,var_dgeQTL))
}
colnames(liste_effects) = c('var_Ad','var_As','corr_Ads','var_Ed','var_Es','corr_Eds','var_C', 'var_dgeQTL')
save(liste_effects, Z, scaled_K, scaled_K_cm, scaled_crossK, scaled_env, scaled_env_cm, scaled_env_cross, scaled_cage_mat, dosages, file = 'data_dir/simulations/polygenicNdgeQTL/image_simulations.RData')


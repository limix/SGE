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



load("output_dir/VD_wE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/all_VCs.RData")
#limit to kind of phenotypes whose correlation Ads we will report, ie those for which aggregate DGE and SGE is high enough (>5%)
all_VCs = all_VCs[all_VCs[,'var_Ad']>0.05 & all_VCs[,'var_As']>0.05,]

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
baseline = c(15,8,0.47,22,16,-0.97,26)
names(baseline) =  c('var_Ad','var_As','corr_Ads','var_Ed','var_Es','corr_Eds','var_C')

soq_var = seq(0,40,by=5)
soq_corr = seq(-1,1,by=0.2)

liste_effects = matrix(NA, nrow = 5 * length(soq_var) + 2 * length(soq_corr), ncol = length(baseline))
colnames(liste_effects) = names(baseline)
roonames = c()
row_counter = 0
for (col in names(baseline)) {
    if (col %in% c('corr_Ads','corr_Eds')) {
        for (j in 1:length(soq_corr)) {
            row_counter = row_counter + 1
            this_baseline = baseline
            this_baseline[col] = soq_corr[j]
            liste_effects[row_counter,] = this_baseline
            roonames = c(roonames,paste(this_baseline,collapse = '_'))
        }
    } else if (col %in% c('var_Ad','var_As','var_Ed','var_Es','var_C')) {
        for (j in 1:length(soq_var)) {
            row_counter = row_counter + 1
            this_baseline = baseline
            this_baseline[col] = soq_var[j]
            liste_effects[row_counter,] = this_baseline
            roonames = c(roonames,paste(this_baseline,collapse = '_'))
        }
    }
}
rownames(liste_effects) = roonames
str(liste_effects)
# num [1:67, 1:7] 0 5 10 15 20 25 30 35 40 15 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : chr [1:67] "0_8_0.47_22_16_-0.97_26" "5_8_0.47_22_16_-0.97_26" "10_8_0.47_22_16_-0.97_26" "15_8_0.47_22_16_-0.97_26" ...
#  ..$ : chr [1:7] "var_Ad" "var_As" "corr_Ads" "var_Ed" ...

save(liste_effects,file='data_dir/simulations/liste_effects.RData')

save.image('data_dir/simulations/image_simulations.RData')


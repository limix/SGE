library(parallel)
load('output_dir/VD_simulations/pruned_dosages_include_DGE_IGE_IEE_cageEffect/all_VCs.RData')
load('data_dir/simulations/liste_effects.RData')

my_f=function(i) {

    load(paste('data_dir/simulations/targets/targets',i,'.RData',sep=''))

    # calculate targets like I did for HS mice
    target_DGE = liste_effects[i,'var_Ad']/total_var * 100
    target_SGE = liste_effects[i,'var_As']/total_var * 100
    target_corr_DGE_SGE = liste_effects[i,'corr_Ads']
    target_DEE = liste_effects[i,'var_Ed']/total_var * 100
    target_SEE = liste_effects[i,'var_Es']/total_var * 100
    target_corr_DEE_SEE = liste_effects[i,'corr_Eds']
    target_cageEffect = liste_effects[i,'var_C']/total_var * 100

    sub = all_VCs[grep(paste('simulation_',rownames(liste_effects)[i],'_',sep=''),all_VCs[,'trait']),]

    res=list(target_DGE = target_DGE,
    DGE = sub[,'var_Ad'],
    target_SGE = target_SGE,
    SGE=sub[,'var_As'],
    target_corr_DGE_SGE = target_corr_DGE_SGE,
    corr_DGE_SGE=sub[,'corr_Ads'],
    target_DEE = target_DEE,
    DEE=sub[,'var_Ed'],
    target_SEE = target_SEE,
    SEE=sub[,'var_Es'],
    target_corr_DEE_SEE = target_corr_DEE_SEE,
    corr_DEE_SEE=sub[,'corr_Eds'],
    target_cageEffect = target_cageEffect,
    cageEffect=sub[,'var_C'],
    target_total_var = target_total_var,
    total_var = sub[,'total_var'])
}

my_try=function(i) {
    tri=try(my_f(i))
    if (inherits(tri,'try-error')) return(NULL) else return(tri)
}
res=mclapply(1:dim(liste_effects)[1],my_try,mc.cores=40)
w=which(unlist(lapply(res,FUN = length))==0)
if (length(w)!=0) stop('pb')


DGE=lapply(res,'[[','DGE')
SGE=lapply(res,'[[','SGE')
corr_DGE_SGE=lapply(res,'[[','corr_DGE_SGE')
DEE=lapply(res,'[[','DEE')
SEE=lapply(res,'[[','SEE')
corr_DEE_SEE=lapply(res,'[[','corr_DEE_SEE')
cageEffect =lapply(res,'[[','cageEffect')
total_var = lapply(res,'[[','total_var')

target_DGE = unlist(lapply(res,'[[','target_DGE'))
target_SGE = unlist(lapply(res,'[[','target_SGE'))
target_corr_DGE_SGE=unlist(lapply(res,'[[','target_corr_DGE_SGE'))
target_DEE=unlist(lapply(res,'[[','target_DEE'))
target_SEE=unlist(lapply(res,'[[','target_SEE'))
target_corr_DEE_SEE=unlist(lapply(res,'[[','target_corr_DEE_SEE'))
target_cageEffect=unlist(lapply(res,'[[','target_cageEffect'))
target_total_var=unlist(lapply(res,'[[','target_total_var'))

pdf('plots_dir/simulations_targetsLikeHSmice.pdf',width=10,height=30)
par(mfrow=c(7,1),mar=c(5,5,1,1))
varying_words = c('Simulated DGE','Simulated SGE','Simulated correlation between DGE and SGE','Simulated DEE','Simulated SEE','Simulated correlation between DEE and SEE','Simulated cage effects')
varying = c('DGE','SGE','corr_DGE_SGE','DEE','SEE','corr_DEE_SEE','cageEffect')
names(varying) = varying_words
width_bp = 0.3

#plotted is value of interest = what will be plotted in boxplots
for (plotted in c('DGE','SGE','corr_DGE_SGE')) {
    
    var_expl = get(plotted)
    #List of 67
    # $ : num [1:500] 0.000369 0.010748 0.003186 0.000119 0.006094 ...
    # $ : num [1:500] 0.0913 0.0858 0.0247 0.0845 0.0569 ...

    target = get(paste('target',plotted,sep='_'))
    # Named num [1:67] 0 3.85 8.2 13.06 18.41 ...
    # - attr(*, "names")= chr [1:67] "var_Ad" "var_Ad" "var_Ad" "var_Ad" ...

    #loop over what is varying
    for (k in 1:length(varying)) {

        vor = varying[k]
        vor_name = names(varying)[k]
        
        # first find range of differences across combinations of effect sizes where where is varying is vor
        if (vor == 'DGE') {from = 1; to = 9}
        if (vor == 'SGE') {from = 10; to = 18}
        if (vor == 'corr_DGE_SGE') {from = 19; to = 29}
        if (vor == 'DEE') {from = 30; to = 38}
        if (vor == 'SEE') {from = 39; to = 47}
        if (vor == 'corr_DEE_SEE') {from = 48; to = 58}
        if (vor == 'cageEffect') {from = 59; to = 67}

        diffs=c()
        for (j in from:to) {
            diffs=c(diffs,var_expl[[j]]-target[[j]])
        }
        length(diffs)
        r=range(diffs,na.rm=T)
        #for legend
        r[2] = r[2] + 0.2*(r[2]-r[1])
        
        ot = 0
        for (j in from:to) {
            #where boxplot will be placed:
            ot = ot + 1

            #put axes labels when first time plotting
            if (j==from) {
                boxplot(var_expl[[j]]-target[[j]],at=j,xlim=c(from-1,to+1),las=1,ylim=r,ylab='',pch=16,cex=0.3,boxwex=width_bp,cex.axis=1.3,cex.lab=1.5)
                mtext(text=paste('Estimated - simulated',plotted),at=mean(r),side=2,padj=-4)
            } else {
                boxplot(var_expl[[j]]-target[[j]],at=j,add=TRUE,pch=16,cex=0.3,boxwex=width_bp,axes=FALSE)
            }
            my_mean=mean(var_expl[[j]]-target[[j]])
            segments(j-0.2,my_mean,j+0.2,my_mean,col='red',lty=1,lwd=3)
        }
        abline(h=0,col='grey')
        mtext(side = 1,text=round(get(paste('target',vor,sep='_'))[from:to], digits = 1), at=from:to,cex=1,padj=2)
        mtext(side = 1,text=vor_name, at=mean(from:to),cex=1,padj=4)
    }
    
}
dev.off()

# baseline
#  var_Ad   var_As corr_Ads   var_Ed   var_Es corr_Eds    var_C 
#   15.00     8.00     0.47    22.00    16.00    -0.97    26.00 
 




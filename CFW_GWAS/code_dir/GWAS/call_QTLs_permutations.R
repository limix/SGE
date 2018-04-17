args <- commandArgs(TRUE)
effect = args[1]
pheno_name = args[2]
alpha = 0.01

library(rhdf5)
task = 'VD_CFW'
kinship_type = 'pruned_dosages'
subset = 'include'

my_get = function(nome) {
	ret = my_h5[[nome]]
	names(ret) = paste(nome,1:length(ret),sep='_')
	return(ret)
}

pvalues_dir = paste('output_dir/pvalues_perms_pointwise_conditional/',effect,'/',kinship_type,'_',subset,'_DGE_IGE_IEE_cageEffect/',sep='')
fole = paste(pvalues_dir,'/',pheno_name,'_','QTLs_alpha',alpha,'.RData',sep='')

if (file.exists(fole)) stop('all done')

h5_file = paste(pvalues_dir,pheno_name,'.h5',sep='')
my_h5 = try(h5read(h5_file,'/'))
if (!inherits(my_h5,'try-error')) {
	my_h5 = as.data.frame(my_h5)
} else my_h5 = NULL

if (is.null(my_h5)) {
	file.remove(h5_file)
	stop('null my_h5')
}

if (any(!c('chr','cumpos','pos') %in% colnames(my_h5))) stop('bad str my_h5')

all_QTLs_perms = NULL
g = grep('pvalues',colnames(my_h5))
if (length(g)==0) stop('no perms')
perms = colnames(my_h5)[g]

nb_perms = length(perms)
for (perm in perms[1:min(100,nb_perms)]) {

	this_h5 = my_h5[,c('chr','pos',perm)]
	this_h5=this_h5[order(this_h5[,perm]),]
	
	if (this_h5[1,perm]>alpha) next

	peaks_pos=c()
	peaks_cumpos=c()
	peaks_chr=c()
	peaks_pvalue=c()
	peaks_ci_starts=c()
	peaks_ci_stops=c()	
	over = F
	while (over==F) {
		qtl_chr=this_h5[1,'chr']
		qtl_pos=this_h5[1,'pos']
		qtl_pvalue=this_h5[1,perm]

		one_side_window1=1500000
	
		qtl_ci_start=qtl_pos-one_side_window1
		qtl_ci_stop=qtl_pos+one_side_window1
		
		peaks_ci_starts=c(peaks_ci_starts,qtl_ci_start)
		peaks_ci_stops=c(peaks_ci_stops,qtl_ci_stop)		
		
		peaks_pos=c(peaks_pos,qtl_pos)
		peaks_chr=c(peaks_chr,qtl_chr)
		peaks_pvalue=c(peaks_pvalue,qtl_pvalue)
		
		remove=which(this_h5[,'chr']==qtl_chr & this_h5[,'pos'] >=qtl_ci_start & this_h5[,'pos']<= qtl_ci_stop)
		this_h5=this_h5[-remove,]		
		this_h5=this_h5[order(this_h5[,perm]),]
		if (dim(this_h5)[1]==0 || this_h5[1,perm]>alpha) over=T
	}
	
	all_peaks=data.frame(perm = perm, measure=pheno_name,chr=peaks_chr,pos=peaks_pos,pvalue = peaks_pvalue, ci_starts=peaks_ci_starts,ci_stops=peaks_ci_stops,stringsAsFactors=F)
	
######
	
	all_peaks$merge=0
	all_new=NULL
	for (chr in unique(all_peaks[,'chr'])) {
		soub=all_peaks[which(all_peaks[,'chr']==chr),]
		changed=T
		while(changed){
			changed=F
			if (dim(soub)[1]!=1) {
				for (k in (1:(dim(soub)[1]-1))) {
					for (l in ((k+1):dim(soub)[1])) {
						if (soub[k,'ci_stops']>=soub[l,'ci_starts'] & soub[k,'ci_starts']<=soub[l,'ci_stops']) {
							changed=T
							if (soub[k,'merge']==0 & soub[l,'merge']==0) {
								soub[k,'merge']=k
								soub[l,'merge']=k
							} else if (soub[k,'merge']!=0 | soub[l,'merge']!=0) {
								soub[k,'merge']=soub[k,'merge']
								soub[l,'merge']=soub[k,'merge']
							}
						}
					}
				}
				
				merge_values=unique(soub[,'merge'])
				if (all(merge_values!=0)) {
					now=NULL
					for (merge_value in merge_values) {
						soubsoub=soub[soub[,'merge']==merge_value,]
						w=which.min(soubsoub[,'pvalue'])
						add=data.frame(perm = perm, measure=pheno_name,chr=chr,pos=soubsoub[w,'pos'],pvalue = soubsoub[w,'pvalue'],ci_starts=min(soubsoub[,'ci_starts']),ci_stops=max(soubsoub[,'ci_stops']),merge=0,stringsAsFactors=F)
						now=rbind(now,add)
					}
					soub=now
				} else {
					now=soub[soub[,'merge']==0,]
					for (merge_value in merge_values[-which(merge_values==0)]) {
						soubsoub=soub[soub[,'merge']==merge_value,]
						w=which.min(soubsoub[,'pvalue'])
						add=data.frame(perm = perm, measure=pheno_name,chr=chr,pos=soubsoub[w,'pos'],pvalue = soubsoub[w,'pvalue'], ci_starts=min(soubsoub[,'ci_starts']),ci_stops=max(soubsoub[,'ci_stops']),merge=0,stringsAsFactors=F)
						now=rbind(now,add)
					}
					soub=now
				}
			}
		}
		all_new=rbind(all_new,soub)
	}
	
	all_QTLs_perms = rbind(all_QTLs_perms,all_new)
}

save(nb_perms,all_QTLs_perms,file = fole)


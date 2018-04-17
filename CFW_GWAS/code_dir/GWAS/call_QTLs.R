library(rhdf5)

pvalues_dir = 'output_dir/pvalues_pointwise_conditional/DGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect'

files = list.files(pvalues_dir)
files = files[grep('h5',files)]
length(files)

pheno2run = read.delim('data_dir/phenotypes/colnames_phenotypes_bc_to_run.txt',as.is = T, header = F)[,1]

inter = intersect(paste(pheno2run,'.h5',sep=''),files)
length(inter)
#170

my_get = function(nome) {
	ret = my_h5[[nome]]
	names(ret) = paste(nome,1:length(ret),sep='_')
	return(ret)
}

alpha = 0.01
all_QTLs = NULL

for (file in inter) {
	pheno_name = sub('.h5','',file,fixed=T)
	print(pheno_name)

	my_h5 = try(h5read(paste(pvalues_dir,'/',pheno_name,'.h5',sep=''),'/'))
	if (inherits(my_h5,'try-error')) next
	if (length(my_h5)!=4 || all(is.na(my_h5$pvalues))) next

	my_h5 = as.data.frame(my_h5)
	my_h5$id = paste(my_h5[,'chr'],my_h5[,'pos'],sep='_')

	my_h5$logP=-log10(my_h5[,'pvalues'])
	my_h5=my_h5[order(my_h5[,'pvalues']),]
	
	if (my_h5[1,'pvalues']>alpha) next

	peaks_marker =c()
	peaks_pos=c()
	peaks_cumpos=c()
	peaks_chr=c()
	peaks_pvalue=c()
	peaks_logP=c()
	peaks_ci_starts=c()
	peaks_ci_stops=c()	
	over = F
	while (over==F) {
		qtl_chr=my_h5[1,'chr']
		qtl_pos=my_h5[1,'pos']
		qtl_cumpos=my_h5[1,'cumpos']
		qtl_marker=my_h5[1,'id']
		qtl_pvalue=my_h5[1,'pvalues']
		qtl_logP=my_h5[1,'logP']
		
		one_side_window1=1500000
	
		qtl_ci_start=qtl_pos-one_side_window1
		qtl_ci_stop=qtl_pos+one_side_window1
				
		peaks_ci_starts=c(peaks_ci_starts,qtl_ci_start)
		peaks_ci_stops=c(peaks_ci_stops,qtl_ci_stop)		
		
		peaks_marker=c(peaks_marker,qtl_marker)
		peaks_pos=c(peaks_pos,qtl_pos)
		peaks_cumpos=c(peaks_cumpos,qtl_cumpos)
		peaks_chr=c(peaks_chr,qtl_chr)
		peaks_pvalue=c(peaks_pvalue,qtl_pvalue)
		peaks_logP=c(peaks_logP,qtl_logP)
		
		remove=which(my_h5[,'chr']==qtl_chr & my_h5[,'pos'] >=qtl_ci_start & my_h5[,'pos']<= qtl_ci_stop)
		my_h5=my_h5[-remove,]		
		my_h5=my_h5[order(my_h5[,'logP'],decreasing=T),]
		if (dim(my_h5)[1]==0 || my_h5[1,'pvalues']>alpha) over=T
	}
	
	all_peaks=data.frame(measure=pheno_name,marker=peaks_marker,chr=peaks_chr,pos=peaks_pos,cumpos=peaks_cumpos,pvalue = peaks_pvalue, logP=peaks_logP,ci_starts=peaks_ci_starts,ci_stops=peaks_ci_stops,stringsAsFactors=F)
	
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
					marker1=soub[k,'marker']
					for (l in ((k+1):dim(soub)[1])) {
						marker2=soub[l,'marker']
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
						w=which.max(soubsoub[,'logP'])
						add=data.frame(measure=pheno_name,marker='merged_peak',chr=chr,pos=soubsoub[w,'pos'],cumpos=soubsoub[w,'cumpos'],pvalue = soubsoub[w,'pvalue'], logP=soubsoub[w,'logP'],ci_starts=min(soubsoub[,'ci_starts']),ci_stops=max(soubsoub[,'ci_stops']),merge=0,stringsAsFactors=F)
						now=rbind(now,add)
					}
					soub=now
				} else {
					now=soub[soub[,'merge']==0,]
					for (merge_value in merge_values[-which(merge_values==0)]) {
						soubsoub=soub[soub[,'merge']==merge_value,]
						w=which.max(soubsoub[,'logP'])
						add=data.frame(measure=pheno_name,marker='merged_peak',chr=chr,pos=soubsoub[w,'pos'],cumpos=soubsoub[w,'cumpos'],pvalue = soubsoub[w,'pvalue'], logP=soubsoub[w,'logP'],ci_starts=min(soubsoub[,'ci_starts']),ci_stops=max(soubsoub[,'ci_stops']),merge=0,stringsAsFactors=F)
						now=rbind(now,add)
					}
					soub=now
				}
			}
		}
		all_new=rbind(all_new,soub)
	}
	all_QTLs=rbind(all_QTLs,all_new)
}

save(all_QTLs,file = paste(pvalues_dir,'/QTLs_alpha',alpha,'.RData',sep=''))




library(rhdf5)

task = 'VD_CFW'
kinship_type = 'pruned_dosages'
subset = 'include'
args <- commandArgs(TRUE)
phenotype=args[1]
locus = args[2]

manhattan_file = paste('plots_dir/fine_maps/fine_maps_',phenotype,'_',locus,'.pdf',sep='')
pdf(manhattan_file, height = 12)
par(mfrow = c(2,1),mar = c(3,6,4,2))
for (effect in c('SGE','DGE')) {
	pvalues_dir = paste('output_dir/fine_map_files/local_pvalues/',effect,'/',kinship_type,'_',subset,'_DGE_IGE_IEE_cageEffect/',sep='')
	file = paste(pvalues_dir,phenotype,'_',locus,'.h5',sep='')
	print(file)
	my_h5 = try(h5read(file,'/'))

	my_h5 = as.data.frame(my_h5,stringsAsFactors = F)
	my_h5$chr = as.numeric(my_h5$chr)
	my_h5$pos = as.numeric(my_h5$pos)
	my_h5$pvalues = as.numeric(my_h5$pvalues)
	my_h5$logP = -log10(my_h5[,'pvalues'])
	my_h5$hwe_fail = as.logical(as.numeric(my_h5$hwe_fail))
	my_h5$info_fail = as.logical(as.numeric(my_h5$info_fail))
	my_h5$pruned_in = as.logical(as.numeric(my_h5$pruned_in))
	my_h5 = my_h5[order(as.numeric(my_h5$pruned_in),my_h5$pos),]

	cols = rep("#053061",dim(my_h5)[1])
	cols[my_h5$info_fail] = 'grey'
	cols[my_h5$hwe_fail] = 'grey'
	cols[my_h5$hwe_fail & my_h5$info_fail] = 'grey'
	cols[my_h5$pruned_in] = "black"

	cexes = rep(0.5,dim(my_h5)[1])
	cexes[my_h5$pruned_in] = 1

	r = range(my_h5[,'pos'])
	chr = unique(my_h5[,'chr'])

	genes=read.delim('data_dir/mart_genes.txt',as.is=T,check.names = F)
	genes=genes[genes[,'Chromosome/scaffold name']==chr & genes[,'Gene start (bp)']<=r[2] & genes[,'Gene end (bp)']>=r[1],]
	genes=genes[order(genes[,'Gene start (bp)']),]
	gms = grepl('^Gm',genes[,'MGI symbol'])
	if (sum(gms)!=0) genes[gms,'MGI symbol'] = 'Gm'
	riks = grepl('Rik$',genes[,'MGI symbol'])
	if (sum(riks)!=0) genes[riks,'MGI symbol'] = 'Rik'
	mirs = grepl('^Mir',genes[,'MGI symbol'])
	if (sum(mirs)!=0) genes[mirs,'MGI symbol'] = 'Mir'
	fams = grepl('^Fam',genes[,'MGI symbol'])
	if (sum(fams)!=0) genes[fams,'MGI symbol'] = 'Fam'
	tmems = grepl('^Tmem',genes[,'MGI symbol'])
	if (sum(tmems)!=0) genes[tmems,'MGI symbol'] = 'Tmem'
	
	range <- max(c(my_h5[,'logP'],5)) + 1
	offset <- ( range * 4 / 3 ) - range
	offset2=0.15

#	plot(my_h5[,c('pos','logP')],pch=16,ylim=c(-offset,range),axes=FALSE,las = 1,ylab = '-logP', cex.lab = 1.5, col =  cols, cex = cexes)	
#	axis(1,cex.axis=1.5,las = 1)
#	axis(2,cex.axis=1.5,ylim = c(0,6),las = 1)
#	title(main=paste(tolower(effect),'GWAS',sep=''),cex.main=1.5)
#	for (i in 1:nrow(genes)) { 
#		k=4
#		if ( genes[i,'Strand'] == 1 ) {
#			arrows(genes[i,'Gene start (bp)'], -offset + (i %% k)*0.5 , genes[i,'Gene end (bp)'], -offset + (i %%  k)*0.5 , length=0.02, lwd=2, code=2, lty="solid",col='red')
#		
#		} else {        
#			arrows(genes[i,'Gene start (bp)'], -offset + (i %%  k)*0.5 , genes[i,'Gene end (bp)'], -offset + (i %%  k)*0.5 , length=0.02, lwd=2, code=1, lty="solid",col='red')
#		}
#		text(genes[i,'Gene start (bp)']+(genes[i,'Gene end (bp)'] - genes[i,'Gene start (bp)'])/2, -offset + (i %% k )*0.5  + 0.2 , labels=genes[i,'MGI symbol'], cex=1.5, las = 2,col='red')
#	}

	#plot without some genes now
	w=which((genes[,'MGI symbol']=='') | gms | riks | mirs | fams | tmems)
	if (length(w)!=0) genes=genes[-w,]
	
	plot(my_h5[,c('pos','logP')],pch=16,ylim=c(-offset,range),axes=FALSE,las = 1,ylab = '-logP', cex.lab = 1.5, col =  cols, cex = cexes)	
	axis(1,cex.axis=1.5,las = 1)
	axis(2,cex.axis=1.5,ylim = c(0,6),las = 1)
	title(main=paste(tolower(effect),'GWAS',sep=''),cex.main=1.5)
	for (i in 1:nrow(genes)) { 
		k=4
		if ( genes[i,'Strand'] == 1 ) {
			arrows(genes[i,'Gene start (bp)'], -offset + (i %% k)*0.5 , genes[i,'Gene end (bp)'], -offset + (i %%  k)*0.5 , length=0.02, lwd=2, code=2, lty="solid",col='red')
		
		} else {        
			arrows(genes[i,'Gene start (bp)'], -offset + (i %%  k)*0.5 , genes[i,'Gene end (bp)'], -offset + (i %%  k)*0.5 , length=0.02, lwd=2, code=1, lty="solid",col='red')
		}
		#do not plot gene names
		#text(genes[i,'Gene start (bp)']+(genes[i,'Gene end (bp)'] - genes[i,'Gene start (bp)'])/2, -offset + (i %% k )*0.5  + 0.2 , labels=genes[i,'MGI symbol'], cex=1.5, las = 2,col='red')
	}

}
dev.off()


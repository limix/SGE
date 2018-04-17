library(RColorBrewer)
library(rhdf5)
task = 'VD_CFW'
kinship_type = 'pruned_dosages'
subset = 'include'
LM = NULL

effect = 'DGE'
pvalues_dir_DGE = paste(paste('output_dir/pvalues_pointwise_conditional',LM,collapse = '_',sep=''),'/',effect,'/',kinship_type,'_',subset,'_DGE_IGE_IEE_cageEffect',sep='')
files_DGE = list.files(pvalues_dir_DGE)
files_DGE = files_DGE[grep('.h5$',files_DGE)]

effect = 'SGE'
pvalues_dir_SGE = paste(paste('output_dir/pvalues_pointwise_conditional',LM,collapse = '_',sep=''),'/',effect,'/',kinship_type,'_',subset,'_DGE_IGE_IEE_cageEffect',sep='')
files_SGE = list.files(pvalues_dir_SGE)
files_SGE = files_SGE[grep('.h5$',files_SGE)]

union = unique(c(files_DGE,files_SGE))
pheno2run = read.delim('data_dir/phenotypes/colnames_phenotypes_bc_to_run.txt',as.is = T, header = F)[,1]
union = intersect(paste(pheno2run,'.h5',sep=''),union)
length(union)
#170

template = h5read(paste(pvalues_dir_SGE,union[1],sep='/'),'/')
template = as.data.frame(template)
template = template[,c(1,3)]

load(paste('output_dir/pvalues_pointwise_conditional/SGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData',sep=''))
SGE_QTLs = all_QTLs
measures2track = unique(SGE_QTLs[,'measure'])
measures2track
# [1] "PAS.First5"   red                "Neuro.DCX_BWcorr"  black          
# [3] "Bioch.Chloride_BWcorr"   blue     "FACS.CD3posCD4pos"    orange       
# [5] "PST.Immobility.First2min"   red  "Hypoxia.f_NR_BWcorr"    purple     
# [7] "Sleep.s12h_D_BWcorr"    red      "BMC.Kurt"          green          
# [9] "Bioch.Calcium_BWcorr"  blue       "PST.Immobility.Last4min"     red
#[11] "Neuro.Ki67_BWcorr"   black         "BMC.osteoporosis"    green        
#[13] "Sleep.VAR_1h_BWcorr"    red      "Sleep.sDif_LD_BWcorr"    red    
#[15] "Cardio.ECG.Tpeak_Tend_BWcorr" brown "Sleep.s12h_L_BWcorr"  red       
#[17] "Haem.abs_neuts"  yellow            
cols_measures2track = c('red','black','blue','orange','red','purple','red','green','blue','red','black','green','red','red','brown','red','yellow')
names(cols_measures2track) = measures2track
sort(cols_measures2track)

load(paste('output_dir/pvalues_pointwise_conditional/DGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData',sep=''))
DGE_QTLs = all_QTLs

draw_manhattan_plot=function(dat,pheno_name,effect,thr=NULL) {
	cols <- c( "#4393C3", "#053061")
	dat$colour=cols[1]
	dat[which(dat[,'chr'] %in% c(2,4,6,8,10,12,14,16,18)),'colour']=cols[2]
	dat$cex = 0.8
	w = which(!(is.na(dat[,'col'])))
	if (length(w)!=0) {
		dat[w,'colour'] = dat[w,'col']
		dat[w,'cex'] = 2
	}
	#otherwise QTL point gets overpainted so is not visible
	dat <- dat[order(dat$logP, decreasing = T),]
	r = range(c(0,range(dat$logP)))
	plot(dat[,c('cumpos','logP')],ylab="",col=dat[,'colour'],xlab="",main='',cex=dat[,'cex'],pch=16,axes=FALSE,xaxs="i", ylim = c(2,10))
	
	title(main=paste('170',effect,'GWAS'),cex.main=2,cex.lab=2)
	axis(2,cex.axis=2,las = 1)

	mtext("-logP",2,padj=-3,cex=2)

	mtext(c(1:19),1,at=tapply(dat$cumpos,dat$chr,median),adj=1,las=1,cex=1.5)
}

all_DGE_h5 = NULL
all_SGE_h5 = NULL
for (pheno_name in sub('.h5','',union)) {
	print(pheno_name)

	DGE_h5 = try(h5read(paste(pvalues_dir_DGE,'/',pheno_name,'.h5',sep=''),'/'))
	if (inherits(DGE_h5,'try-error') | length(DGE_h5)!=4) next
	if (all(is.na(DGE_h5$pvalues))) stop('pb')
	DGE_h5 = as.data.frame(DGE_h5,stringsAsFactors =F)
	if (any(DGE_h5$chr != template[,1] | DGE_h5$pos != template[,2])) stop('pb markers positions')
	DGE_h5$chr = as.numeric(sub('chr','',DGE_h5$chr))
	DGE_h5$logP=-log10(DGE_h5$pvalues)
	#DGE_h5=DGE_h5[order(DGE_h5[,'logP'],decreasing=T),]
	H5close()
	DGE_h5$col = NA
	w = which(DGE_QTLs[,'measure'] == pheno_name)
	for (i in w) {
		quels = which(DGE_h5[,'chr'] == DGE_QTLs[i,'chr'] & DGE_h5[,'pos'] == DGE_QTLs[i,'pos'])
		if (length(quels)== 0) stop() 
		if (pheno_name %in% measures2track) DGE_h5[quels,'col'] = cols_measures2track[pheno_name] else DGE_h5[quels,'col'] = 'darkgrey'
	}
	all_DGE_h5 = rbind(all_DGE_h5, DGE_h5[DGE_h5[,'logP']>=2,])

	SGE_h5 = try(h5read(paste(pvalues_dir_SGE,'/',pheno_name,'.h5',sep=''),'/'))
	if (inherits(SGE_h5,'try-error')| length(SGE_h5)!=4) next
	if (all(is.na(SGE_h5$pvalues))) stop('pb')
	if (any(SGE_h5$chr != template[,1] | SGE_h5$pos != template[,2])) stop('pb markers positions')
	SGE_h5 = as.data.frame(SGE_h5,stringsAsFactors =F)
	SGE_h5$chr = as.numeric(sub('chr','',SGE_h5$chr))
	SGE_h5[,'logP']=-log10(SGE_h5[,'pvalues'])
	#SGE_h5=SGE_h5[order(SGE_h5[,'logP'],decreasing=T),]
	H5close()
	SGE_h5$col = NA
	w = which(SGE_QTLs[,'measure'] == pheno_name)
	if (length(w) == 0) DGE_h5[,'in_QTLs'] = NA
	for (i in w) {
		quels = which(SGE_h5[,'chr'] == SGE_QTLs[i,'chr'] & SGE_h5[,'pos'] == SGE_QTLs[i,'pos'])
		if (length(quels)== 0) stop() 
		if (pheno_name %in% measures2track) SGE_h5[quels,'col'] = cols_measures2track[pheno_name] else SGE_h5[quels,'col'] = 'darkgrey'
	}
	all_SGE_h5 = rbind(all_SGE_h5, SGE_h5[SGE_h5[,'logP']>=2,])

}

all_DGE_h5[all_DGE_h5[,'logP']>10,'logP'] = 10
w = which(all_DGE_h5[,'col'] == 'black')
if (length(w)!=0) all_DGE_h5[w,'col'] = 'darkgrey'

jpeg('plots_dir/manhattans_pointwise_conditional/porcupine.jpg',width=1300,height = 1000)
par(mfrow=c(2,1),mar = c(3,6,4,2))
draw_manhattan_plot(all_SGE_h5,pheno_name,effect = 'SGE')
draw_manhattan_plot(all_DGE_h5,pheno_name,effect = 'DGE')
graphics.off()


effect = 'DGE'
load(paste('output_dir/pvalues_pointwise_conditional/',effect,'/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData',sep=''))
all_QTLs = all_QTLs[,c('measure','chr','pos','logP','FDR')]
colnames(all_QTLs) = c('Name','Chr','Pos (bp)','logP','Per-phenotype Q value')
all_QTLs[,'logP'] = round(all_QTLs[,'logP'],digits = 1)
all_QTLs[,'Per-phenotype Q value'] = round(all_QTLs[,'Per-phenotype Q value'],digits = 2)

load('output_dir/effect_sizes/betas_conditioning/var_expl/betas.RData')
motch = match(paste(all_QTLs[,'Name'],'_','chr',all_QTLs[,'Chr'],'_',all_QTLs[,'Pos (bp)'],'_','direct',sep=''),paste(indiv_DGE_QTLs_betas[,1],indiv_DGE_QTLs_betas[,2],sep='_'))
any(is.na(motch))
#FALSE

all_QTLs$effect_size = round(indiv_DGE_QTLs_betas[motch,3],digits = 1)
colnames(all_QTLs) = sub('effect_size','% variance',colnames(all_QTLs))

model_menu = read.delim('data_dir/phenotypes/TableS1_phenotypes.txt',as.is=T,check.name = F)
model_menu = model_menu[1:200,]

motch = match(sub('_BWcorr','',all_QTLs[,'Name']), model_menu[,'Name'])
any(is.na(motch))
#TRUE - osteoporosis manually changed later

all_QTLs$Phenotype = model_menu[motch,'Measure']

all_QTLs = all_QTLs[,c('Name','Chr','Pos (bp)','logP','% variance')]

all_QTLs = all_QTLs[order(all_QTLs[,'Chr'],all_QTLs[,'Pos (bp)']),]
write.table(all_QTLs,file = 'plots_dir/DGE_QTLs.xls',sep='\t',quote = F, row.names = F, col.names = T)


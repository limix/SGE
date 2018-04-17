library(ff)

args <- commandArgs(TRUE)
chr = args[1]
beg = args[2]
end = args[3]

load(paste("data_dir/old02202018_dosages/STITCH/EMparameters.chr",chr,".RData",sep=""))
#!!! that loads chr = "chr16" (when it used to be 16), also loads a "K" and "T"
dosages.ff=ff(filename=paste("data_dir/old02202018_dosages/STITCH/dosages.",chr,".ff",sep=""),vmode="double",dim=c(T,N),readonly=TRUE)
rm(T,K)
#these are between 0 and 1 - multiply by 2 (Robbies Davies, personal communication)
dosages = as.matrix(as.data.frame(as.ffdf(dosages.ff)))*2

load(paste("data_dir/old02202018_dosages/STITCH/sampleNames.",chr,".RData",sep=""))
sampleNames = sub('_recal.*','',sampleNames,perl=T)
sampleNames = sub('___','/',sampleNames,fixed=T)

load(paste('data_dir/old02202018_dosages/STITCH/pos.',chr,'.RData',sep=''))
map = data.frame(sub('chr','',as.character(pos[,'CHR'])),pos[,'POS'],pos[,'POS'],paste(as.character(pos[,'REF']),as.character(pos[,'ALT']),sep='/'),'+',paste(as.character(pos[,'CHR']),as.character(pos[,'POS']),sep='_'),stringsAsFactors = F)
colnames(map) = c('chr','pos','pos2','allele','strand','id')
load(paste('data_dir/old02202018_dosages/STITCH/info.',chr,'.RData',sep=''))
load(paste('data_dir/old02202018_dosages/STITCH/hwe.',chr,'.RData',sep=''))

w = which((map[,'pos']>=as.numeric(beg)) & (map[,'pos']<=as.numeric(end)))
length(w)

dosages = dosages[w,]
map = map[w,]
hwe = hwe[w]
info = info[w]

hwe_fail = as.integer(hwe<1e-6)
info_fail = as.integer(info<0.4)

pruned = read.table('data_dir/dosages/pruned_dosages/pruned_SNPs_IDs.txt',as.is=T)[,1]
pruned_in = as.integer(map[,'id'] %in% pruned)

dir.create('output_dir/fine_map_files/all_local_dosages/',recursive = T)
save(dosages,map,hwe_fail,info_fail,pruned_in,sampleNames,file = paste('output_dir/fine_map_files/all_local_dosages/',chr,'_',beg,'_',end,'.RData',sep=''))


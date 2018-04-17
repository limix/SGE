load('data_dir/dosages/pruned_dosages/std_my_final_dosages.RData')
dim(std_dosages)
#[1]   1869 353697

library(rhdf5)
fid='data_dir/CFWmice.h5'
h5createGroup(fid,"direct_pruned_dosages")
h5createGroup(fid,"direct_pruned_dosages/row_header")
max(nchar(rownames(std_dosages)))
#15
h5createDataset(file=fid,dataset="direct_pruned_dosages/row_header/sample_ID",dims=dim(std_dosages)[1],storage.mode='character',size=20)
h5write(obj=rownames(std_dosages),file=fid,name="direct_pruned_dosages/row_header/sample_ID")

pos = do.call('rbind',sapply(colnames(std_dosages),strsplit,'_'))
pos = as.data.frame(pos,stringsAsFactors=F)
pos[,1] = as.numeric(sub('chr','',pos[,1]))
pos[,2] = as.numeric(pos[,2])
pos = as.matrix(pos)
ordre = order(pos[,1],pos[,2])
all(ordre == seq(1,dim(pos)[1]))
#TRUE
cumpos = pos[pos[,1]==1,2]
for (chr in unique(pos[,1])[-1]) {
	w = which(pos[,1]==chr)
	cumpos = c(cumpos,pos[w,2]+cumpos[length(cumpos)])
}
#plot(pos[,1],cumpos)
ordre = order(cumpos)
all(ordre == seq(1,length(cumpos)))
#TRUE
h5createGroup(fid,"direct_pruned_dosages/col_header")
h5write(obj=pos[,1],file=fid,name="direct_pruned_dosages/col_header/chr")
h5write(obj=pos[,2],file=fid,name="direct_pruned_dosages/col_header/pos")
h5write(obj=cumpos,file=fid,name="direct_pruned_dosages/col_header/cumpos")
h5write(obj=colnames(std_dosages),file=fid,name="direct_pruned_dosages/col_header/id")

h5createDataset(file=fid,dataset="direct_pruned_dosages/matrix",dims=dim(std_dosages),chunk = c(dim(std_dosages)[1],ceiling(dim(std_dosages)[2])/2))
h5write(obj=std_dosages,file=fid,name="direct_pruned_dosages/matrix")


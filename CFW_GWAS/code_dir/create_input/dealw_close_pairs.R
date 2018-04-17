load('data_dir/dosages/pruned_dosages/GRM.RData')

library(MASS)
truehist(GRM[lower.tri(GRM,diag=F)],ylim=c(0,0.01),proba=TRUE, nbins = 10)
#heavy right tail. 
#based on distribution, will use 0.3 as cutoff to define close pairs
close_pairs=NULL
for (i in 1:dim(GRM)[1]) {
	w=which(GRM[i,]>0.3)
	w=w[w>i]
	if (length(w)!=0) {
		for (k in w) {
			close_pairs = rbind(close_pairs,c(colnames(GRM)[i],colnames(GRM)[k],GRM[i,k]))
		}
	}
}
dim(close_pairs)
length(unique(c(close_pairs[,1],close_pairs[,2])))
#  mice involved in close pairs thus defined

#bring in cage information
load('data_dir/phenotypes/imputed_cage_info.RData')
cages = imputed_cage_info[,1]
length(unique(cages))
#[1] 716
names(cages) = rownames(imputed_cage_info)
close_pairs = as.data.frame(close_pairs,stringsAsFactors = F)
motch = match(close_pairs[,1],names(cages))
any(is.na(motch))
#FALSE
close_pairs$cage1 = cages[motch]
motch = match(close_pairs[,2],names(cages))
any(is.na(motch))
#FALSE
close_pairs$cage2 = cages[motch]

length(unique(c(close_pairs[,'cage1'],close_pairs[,'cage2'])))
#[1] 145 out of 716 cages involved to some extent in close pairs (look at cages rather than mice here because will exclude all mice in a cage where one of mice excluded)

table(table(c(close_pairs[,'cage1'],close_pairs[,'cage2'])))
# 1  2  3  4  5  6  9 
#83 26 17 10  5  3  1 
# so some cages show up in many close pairs
sum(table(c(close_pairs[,'cage1'],close_pairs[,'cage2']))>3)
#19

toble = table(c(close_pairs[,'cage1'],close_pairs[,'cage2']))
pbic_cages = names(toble)[toble>3]
length(pbic_cages)
#19
all_cagemates_close_pairs = names(cages)[cages %in% pbic_cages]
length(all_cagemates_close_pairs)
#57
save(all_cagemates_close_pairs,file = 'data_dir/dosages/close_pairs_GRM.RData')

library(rhdf5)
fid='data_dir/CFWmice.h5'

motch = match(all_cagemates_close_pairs, rownames(GRM))
any(is.na(motch))
#FALSE
include = rownames(GRM)[-motch]
length(include)
#[1] 1812

h5createGroup(file=fid,'subsets')
max(nchar(include))
#15
h5createDataset(file=fid,dataset="subsets/include",dims=length(include),storage.mode='character',size=20)
h5write(obj=include,file=fid,name="subsets/include")




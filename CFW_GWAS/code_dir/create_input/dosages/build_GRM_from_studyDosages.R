library(parallel)

load('data_dir/dosages/pruned_dosages/std_my_final_dosages.RData')

tot = dim(std_dosages)[2]
by = 2000
mox = ceiling(tot/by)
my_test = function(start) {
	if (start==mox) r = seq((start-1)*by+1,tot) else r = seq((start-1)*by+1,(start-1)*by+by)
	return(r)
}
test = unlist(mclapply(1:mox,my_test,mc.cores=12))
all((1:tot) == sort(test))
#[1] TRUE

tot = dim(std_dosages)[2]
by = 2000
mox = ceiling(tot/by)
my_do = function(start) {
	print(start)
	if (start==mox) r = seq((start-1)*by+1,tot) else r = seq((start-1)*by+1,(start-1)*by+by)
	GRM = std_dosages[,r] %*% t(std_dosages[,r])
	return(GRM)
}
GRMs = mclapply(1:mox,my_do,mc.cores=12)
for (seg in 1:mox) {
	if (seg==1) GRM = GRMs[[1]] else GRM = GRM + GRMs[[seg]]
}
dim(GRM)
GRM = GRM/dim(std_dosages)[2]
#as long as dosages stdized, dividing by number of SNPs is equivalent to dividing by sum of diagonal

rownames(GRM) = colnames(GRM) = rownames(std_dosages)

save(GRM, file = 'data_dir/dosages/pruned_dosages/GRM.RData')


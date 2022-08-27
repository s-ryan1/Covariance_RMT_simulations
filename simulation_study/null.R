library(cpt.cov)
library(purrr)

#args = commandArgs(trailingOnly=TRUE)
#iter = as.numeric(args[1])

run.sim = function(i, n, p){
	set.seed(353*n + 541*p + 7*i)
    data = rnorm(n*p,0,1) %>% matrix(ncol=p) 		
    return(max(matrix.dist.test.stat(data, 4*p), na.rm=T))
}

n = c(200,500,1000,2000,5000)
p= c(10,20,50,100)
run = seq(0,9)*100

#p = 100

parameters = cross3(n,p,run) 
n= parameters[[iter]][[1]]
p= parameters[[iter]][[2]]
run= parameters[[iter]][[3]]
result=map(seq(1,100), ~run.sim(run+.x, n, p))


string = sprintf("~/Documents/covariance/null/Results/null%d_%d_%d.Rds",  n, p, run) 
saveRDS(result, string)

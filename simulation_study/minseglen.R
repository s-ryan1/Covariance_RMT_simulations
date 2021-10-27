library(cpt.cov)
library(tibble)
library(tidyr)
library(purrr)

args = commandArgs(trailingOnly=TRUE)
iter = as.numeric(args[1])

run.sim = function(i, n, p){
	set.seed(353*n + 541*p + 7*i)
    data = rnorm(n*p,0,1) %>% matrix(ncol=p) 		
    out=matrix.dist.test.stat(data,1)
    minseglens = c(1.1*p,1.2*p,1.5*p,2*p,4*p,8*p,30)
    minseglens = floor(minseglens)
    values = map_dbl(minseglens, ~max(out[(.x+1):(n-.x)],na.rm=T))
    out = tibble(n=n,p=p,minseglens=minseglens,values=values) 
    return(out)
}

params = list(c(100,3), c(500,15), c(2500,100))
run = seq(0,9)*100

#p = 100

n=params[[iter%%3+1]][1]
p=params[[iter%%3+1]][2]
run=run[iter%%10+1]
result=map(seq(1,100), ~run.sim(run+.x, n, p)) %>% tibble %>% unnest


string = sprintf("~/Documents/covariance/minseglen/Results/updated%d.Rds", iter) 
saveRDS(result, string)

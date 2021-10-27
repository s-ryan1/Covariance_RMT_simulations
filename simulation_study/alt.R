library(cpt.cov)
library(purrr)

#args = commandArgs(trailingOnly=TRUE)
#iter = as.numeric(args[1])
run.sim = function(i, n, p, delta){

    #scale = list(diag(p), delta*diag(p))
    #rotation = list(diag(p), delta*B) 

	set.seed(353*n + 541*p + 7*i)
    noise = rnorm(n*p,0, 1) %>% matrix(ncol=p)
    noise = list(noise[1:(n/2),], delta*noise[-c(1:n/2),])
    #noise = map2(noise, wish, ~.x%*%sqrtm(.y))

    data = noise %>% reduce(rbind)
    out = matrix.dist.test.stat(data, 4*p)
    return(c(max(out, na.rm=T), which.max(out)))
}



n = c(200,500,1000,2000,5000)
p= c(10,50,100)
run = seq(0,9)*100
delta=c(1.05,1.1,1.15,1.2)

#p = 100

parameters = cross(list(n,p,delta,run)) 
n= parameters[[iter]][[1]]
p= parameters[[iter]][[2]]
delta=parameters[[iter]][[3]]
run= parameters[[iter]][[4]]
result=map(seq(1,100), ~run.sim(run+.x, n, p,delta))

#result = list(n=n,p=p,result=result) 





string = sprintf("~/Documents/covariance/alt/Results/alt%d.Rds", iter) 
saveRDS(result, string)

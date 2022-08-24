library(cpt.cov)


bonferoni = function(n,alpha){ qnorm(1-.05/n)}

n=500
p=30
Y = matrix(rnorm(n*p,0,1), nrow=n)
Y[300:500,] = 2*Y[300:500,]

minseglen = 120
result =  matrix.dist.test.stat(Y, minseglen-p)

if(max(result, na.rm=T)>bonferoni(n,.05)){
    print(sprintf("Change detected at %d", which.max(abs(result))))
} else {
    print(sprintf("No change detected"))
} 

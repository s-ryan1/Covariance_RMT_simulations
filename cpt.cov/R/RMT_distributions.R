marchenko.pastur = function(lambda){ 
    a=(1-sqrt(lambda))^2
    b=(1+sqrt(lambda))^2
    func = function(x){
       return((log(x)^2)*(1/(2*pi*lambda*x))*sqrt((x-a)*(b-x)))
    }
    return(list(func,a,b))
}

#' @export
fisher.esd = function(y1,y2){
    c2=y2#max(c(y1,y2))
    c1=y1#min(c(y1,y2))
    h = sqrt(c1+c2-c1*c2)
    a = ((1-h)^2)/(1-c2)^2
    b = ((1+h)^2)/(1-c2)^2
    function(x){
        result = rep(0, length(x))
        result[x>a & x<b] = ((1-c2)*sqrt((b-x)*(x-a)))/(2*pi*x*(c1 +x*c2))
        return(result)
    }
}

#' @export
construct.fisher.support= function(y1,y2){
    c2=y2#max(c(y1,y2))
    c1=y1#min(c(y1,y2))
    h=sqrt(c1+c2-c1*c2)
    a = ((1-h)^2)/(1-c2)^2
    b = ((1+h)^2)/(1-c2)^2
    return(c(a,b))
}

#' @export
fisher.integrator=function(f,n,p,epsilon){
    function(t){
            integrand = function.prod(f, fisher.esd(p/t, p/(n-t)))
            support = construct.fisher.support(p/t, p/(n-t))
            result =integrate(integrand, support[[1]], support[[2]])[[1]]
            return(result)
    }
}

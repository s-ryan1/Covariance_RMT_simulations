#' @export
bin.seg = function(data, parent, method, threshold, minseglen, params){

    if(nrow(data) > 2*(minseglen+ncol(data))){
        test.stat = exec(method, data, minseglen, !!!params)
        candidate = which.max(test.stat)
        if(test.stat[candidate] > threshold){
            data = list(data[seq(1, candidate),], data[seq(candidate+1, nrow(data)),])
            parents = list(c(parent[1],parent[1] + candidate), c(parent[1] +candidate, parent[2]))
            cpt = map2(data, parents, bin.seg, method, threshold, minseglen, params)
            #cpt = reduce(cpt, rbind)
            return(list(c(parent, test.stat[candidate], candidate+parent[1]), cpt[[1]], cpt[[2]]))
        } else {
            return(list(c(parent, test.stat[candidate], candidate+parent[1]), 0))
        }
    }else{
        return(list(c(parent, NA, NA), 0))
    }

}

#' @export
bin.seg.to.cpt = function(result, threshold){
    cpts = c()
    cpts.out = c()
    thresholds = c()
    a = function(segment){
        if(!is.na(segment[[1]][3])){
           if(segment[[1]][3] > threshold){
               cpts <<- rbind(cpts, segment[[1]]) 
               thresholds <<- c(thresholds,segment[[1]][3])
               walk(segment[-1], a)
           } else {
               cpts.out <<- rbind(cpts.out, segment[[1]]) 
           }
        } else {
           cpts.out <<- rbind(cpts.out, segment[[1]]) 
        }
    }
    a(result)
    
    cpts.out[,1] = cpts.out[,1]+1

    return(list(cpts=cpts[,4],segments = cpts.out, thresholds = thresholds))
}

#' @export
minseglen.max = function(X, minseglen){
    return(max(X[seq(floor(minseglen)+1,length(X)-floor(minseglen))], na.rm=TRUE))
}

#' @export
detection.rates = function(cpts.est, cpts.true, detect.interval){
    cpts = cross2(cpts.est, cpts.true)
    cpts = discard(cpts, ~abs(.x[[1]] - .x[[2]]) > detect.interval)
    correct.locs = unique(map_dbl(cpts, ~.x[[2]]))
    tdr = 0; fdr=0;
    if(length(cpts.true)>0 & length(correct.locs)>0){
            tdr = length(correct.locs)/length(cpts.true)
    }
    if(length(cpts.est)>0 & length(cpts.true)>0){
            fdr = (length(cpts.est)-length(correct.locs))/length(cpts.est)
    if(fdr < 0){ print(list(cpts.est,cpts.true))}
    }
    #if(length(cpts.est)==0){ 
    #    correct.locs=NA
    #    if(length(cpts.true) == 0 ){ 
    #        tdr = 1; fdr=0
    #    } else if(length(cpts.true)>0) {
    #        tdr=0; fdr=0
    #    }
    #}

    return(list(m= abs(length(cpts.est) - length(cpts.true)), correct.locs = correct.locs, TDR=tdr, FDR=fdr))
}

#' @export
MAE = function(cpts.est, cpts.true, data){
    segments = c(0,unique(sort(c(cpts.est, cpts.true))), nrow(data))
    segments = imap(segments[-1], ~c(segments[.y]+1, segments[.y+1])) 
    return(sum((map_dbl(segments, interval.mae, cpts.est, cpts.true, data)))/nrow(data))
}

#' @export 
spectral.error = function(cpts.est, cpts.true, data){
    segments = c(0,unique(sort(c(cpts.est, cpts.true))), nrow(data))
    segments = imap(segments[-1], ~c(segments[.y]+1, segments[.y+1])) 
    return(sum((map_dbl(segments, interval.spectral.error, cpts.est, cpts.true, data)))/nrow(data))
}

interval.mae = function(interval, cpts.est, cpts.true, data){
    n=nrow(data)
    cpts.true = c(0,cpts.true, n)
    cpts.est = c(0,cpts.est, n)
    left.true = max(cpts.true[cpts.true<interval[1]])
    left.est = max(cpts.est[cpts.est<interval[1]])
    right.true = min(cpts.true[cpts.true>=interval[2]])
    right.est = min(cpts.est[cpts.est>=interval[2]])
    return((diff(interval)+1)*(sum(abs(cov(data[left.est:right.est,]) - cov(data[left.true:right.true,])))))
}

interval.spectral.error = function(interval, cpts.est, cpts.true, data){
    n=nrow(data)
    cpts.true = c(0,cpts.true, n)
    cpts.est = c(0,cpts.est, n)
    left.true = max(cpts.true[cpts.true<interval[1]])
    left.est = max(cpts.est[cpts.est<interval[1]])
    right.true = min(cpts.true[cpts.true>=interval[2]])
    right.est = min(cpts.est[cpts.est>=interval[2]])
    return((diff(interval)+1)*(sum(abs(eigen(cov(data[left.est:right.est,]))[[1]] - eigen(cov(data[left.true:right.true,]))[[1]]))))
}

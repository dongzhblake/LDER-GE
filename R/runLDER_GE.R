#' @title Main function
#' @description Run LDER-GE 
#' @param assoc GWIS summary statistics, need to include snp, chr, a0, a1, z (header is necessary)
#' @param n.gwas The sample size of the GWAS summary statistics
#' @param path The path of LD panel directory
#' @param LD.insample T/F, whether the LD reference is estimated with the GWAS cohort (T) or external reference panel (e.g. 1000 Genome Project) (F)
#' @param n.ld The sample size of the LD reference
#' @param cores The number of cores for computation in parallel
#' @param method 'lder', 'ldsc', or 'both'
#' @param size_num Number of blocks for jackknife
#' @param type if 'jack' then conducts delete-wise block jackknife. if 'none' then does not do inference 
#' @import  data.table stats utils
#' @export
#'
#'
runLDER_GE <- function(assoc, n.gwas, path, LD.insample=T,  n.ld,method='lder', cores=10, a=NULL, type='jack',size_num=200){
  A=unlist(strsplit(list.files(path,pattern="INFO"),".txt"));B=unlist(strsplit(A,"SNPINFO"));C=max(as.numeric(B),na.rm=T)
  library(parallel)
  ldpath <- ldpath.shrink <- path
  print("matching summary statistics with LD panel")
  stats <- mclapply(0:C,get.stats,assoc=assoc,ldpath=ldpath,ldpath.shrink=ldpath.shrink,mc.cores=cores,n.ld=n.ld)
  if(method=='lder'){
    res <- lder(stats=stats,n.gwas=n.gwas,a=NULL,rough=!LD.insample,cores=cores,twostage=!LD.insample,type=type)
    return(res)
  }else if(method=='ldsc'){
    res <- ldsc(stats=stats,n.gwas=n.gwas,a=NULL,cores=cores,twostage=F,type=type)
    return(res)
  }else if(method=='both'){
    res1 <- lder(stats=stats,n.gwas=n.gwas,a=NULL,rough=!LD.insample,cores=cores,twostage=!LD.insample,type=type)
    res2 <- ldsc(stats=stats,n.gwas=n.gwas,a=NULL,cores=cores,twostage=F,type=type)
    return(list(lder=res1,ldsc=res2))
  }
}

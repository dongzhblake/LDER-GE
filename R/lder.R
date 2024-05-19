
lder <- function(stats,n.gwas,a=NULL,rough=F,cores=20,twostage=T,type='jack',size_num=200){
  library(parallel)
  temp <- unlist(lapply(stats,length))
  stats[temp==1] <- NULL
  x1 <- unlist(sapply(stats,'[[',1))
  lam1 <- unlist(sapply(stats,'[[',3))
  do.jack <- function(stats,j,cores,n.gwas,a,rough,twostage){
    stats1 <- stats[-j]
    x1 <- unlist(sapply(stats1,'[[',1))
    lam1 <- unlist(sapply(stats1,'[[',3))
    res <- get.res(x1,lam1,n.gwas=n.gwas,a=a,rough=rough,twostage=twostage)
    res$h2 <- res$h2/length(x1)
    return(res)
  }
  res <- get.res(x1,lam1,n.gwas=n.gwas,a=a,rough=rough,twostage=twostage)
  if(type=='none'){
    return(list(h2I=res$h2,intecept=res$a*n.gwas+1 ))
  }
  
  mj=sapply(1:length(stats),FUN=function(x)(length(stats[[x]]$ldsc)))
  sizedf=as.data.frame(cbind(mj,1:length(stats)))
  define_block_unit <- function(size_num=200){
    sizes_perb=c()
    for(b in seq(0.5,2.5,by=0.01)){
      target_size=sum(mj)/size_num/b
      grouped_block=list()
      size_each=c()
      size_now=0
      block_group=c()
      for(i in 1:nrow(sizedf)){
        if(size_now<target_size){
          size_now=size_now+sizedf$mj[i]
          block_group=c(block_group,sizedf$V2[i])}
        else{
          grouped_block=append(grouped_block,list(block_group))
          size_each=c(size_each,size_now)
          size_now=sizedf$mj[i]
          block_group=c(sizedf$V2[i])}}
      sizes_perb=c(sizes_perb,length(grouped_block))}
    b=seq(0.5,2.5,by=0.01)[which(abs(sizes_perb-size_num)==min(abs(sizes_perb-size_num)))][1]
    target_size=sum(mj)/size_num/b
    grouped_block=list()
    size_each=c()
    size_now=0
    block_group=c()
    for(i in 1:nrow(sizedf)){
      if(size_now<target_size){
        size_now=size_now+sizedf$mj[i]
        block_group=c(block_group,sizedf$V2[i])}
      else{
        grouped_block=append(grouped_block,list(block_group))
        size_each=c(size_each,size_now)
        size_now=sizedf$mj[i]
        block_group=c(sizedf$V2[i])}}
    sizes_perb=c(sizes_perb,length(grouped_block))
    return(return(grouped_block))
  }
  grouped_block=define_block_unit(size_num=size_num)
  
  if(type=='jack'){
    pb <- txtProgressBar(0,length(grouped_block),style=3)
    temp1.jack=list()
    print("Estimating SE with delete-block-jackknife")
    for(t in 1:length(grouped_block)){temp1.jack=append(temp1.jack,list(do.jack(stats=stats,j=grouped_block[[t]],cores=cores,n.gwas=n.gwas,a=a,rough=rough,twostage=F)))
    setTxtProgressBar(pb, t)}
    print("delete-block-jackknife done")
    close(pb)
    hh2 <-  unlist(sapply(temp1.jack,'[[','h2'))
    inteceptf2 <- unlist(sapply(temp1.jack,'[[','a'))
    n.block <- length(grouped_block)
    h33.jack <- sd(hh2)*sqrt(n.block)*length(x1)
    a33.jack <- sd(inteceptf2)*sqrt(n.block)
    return(list(h2I=res$h2,intecept=res$a*n.gwas+1,h2I.se=h33.jack,h2I.p=pchisq((res$h2/h33.jack)^2,df=1,lower.tail = F),intecept.se=a33.jack*n.gwas))
  }
}


ldsc <- function(stats,n.gwas,a=NULL,cores=20,twostage=T,type='jack',size_num=200){
  library(parallel)
  temp <- unlist(lapply(stats,length))
  stats[temp==1] <- NULL
  z1 <- unlist(lapply(stats,'[[','z'))
  ldsc1 <- unlist(sapply(stats,'[[','ldsc'))
  do.jack.ldsc <- function(stats,j,cores,n.gwas,a,rough,twostage){
    stats1 <- stats[-j]
    z1 <- unlist(sapply(stats1,'[[','z'))
    ldsc1 <- unlist(sapply(stats1,'[[','ldsc'))
    res <- get.res.ldsc(z1,ldsc1,n.gwas=n.gwas,a=a,twostage=twostage)
    res$h2 <- res$h2/length(z1)
    return(res)
  }
  res <- get.res.ldsc(z1,ldsc1,n.gwas=n.gwas,a=a,twostage=twostage)
  if(type=='none'){
    return(list(h2I=res$h2,intecept=res$a*n.gwas+1))
  }
  
  mj=sapply(1:length(stats),FUN=function(x)(length(stats[[x]]$ldsc)))
  sizedf=as.data.frame(cbind(mj,1:length(stats)))
  define_block_unit <- function(size_num=200){
    sizes_perb=c()
    for(b in seq(0.5,2.5,by=0.01)){
      target_size=sum(mj)/size_num/b
      grouped_block=list()
      size_each=c()
      size_now=0
      block_group=c()
      for(i in 1:nrow(sizedf)){
        if(size_now<target_size){
          size_now=size_now+sizedf$mj[i]
          block_group=c(block_group,sizedf$V2[i])}
        else{
          grouped_block=append(grouped_block,list(block_group))
          size_each=c(size_each,size_now)
          size_now=sizedf$mj[i]
          block_group=c(sizedf$V2[i])}}
      sizes_perb=c(sizes_perb,length(grouped_block))}
    b=seq(0.5,2.5,by=0.01)[which(abs(sizes_perb-size_num)==min(abs(sizes_perb-size_num)))][1]
    target_size=sum(mj)/size_num/b
    grouped_block=list()
    size_each=c()
    size_now=0
    block_group=c()
    for(i in 1:nrow(sizedf)){
      if(size_now<target_size){
        size_now=size_now+sizedf$mj[i]
        block_group=c(block_group,sizedf$V2[i])}
      else{
        grouped_block=append(grouped_block,list(block_group))
        size_each=c(size_each,size_now)
        size_now=sizedf$mj[i]
        block_group=c(sizedf$V2[i])}}
    sizes_perb=c(sizes_perb,length(grouped_block))
    return(return(grouped_block))
  }
  grouped_block=define_block_unit(size_num=size_num)
  
  
  if(type=='jack'){
    pb <- txtProgressBar(0,length(grouped_block),style=3)
    temp1.jack=list()
    print("Estimating SE with delete-block-jackknife")
    for(t in 1:length(grouped_block)){temp1.jack=append(temp1.jack,list(do.jack.ldsc(stats=stats,j=grouped_block[[t]],cores=cores,n.gwas=n.gwas,a=a,rough=rough,twostage=F)))
    setTxtProgressBar(pb, t)}
    print("delete-block-jackknife done")
    close(pb)
    hh2 <-  unlist(sapply(temp1.jack,'[[','h2'))
    inteceptf2 <- unlist(sapply(temp1.jack,'[[','a'))
    n.block <- length(grouped_block)
    h33.jack <- sd(hh2)*sqrt(n.block)*length(z1)
    a33.jack <- sd(inteceptf2)*sqrt(n.block)
    return(list(h2I=res$h2,intecept=res$a*n.gwas+1,h2I.se=h33.jack,h2I.p=pchisq((res$h2/h33.jack)^2,df=1,lower.tail = F),intecept.se=a33.jack*n.gwas))
  }
}




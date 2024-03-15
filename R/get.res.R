get.res <- function(x1,lam1,n.gwas,a,rough,twostage){
  m <- length(x1)
  if(twostage){
    s2thld <- quantile(x1^2, 0.9909)
    idx <- (x1^2<=s2thld)
    temp1 <-  calH2.new1(x1[idx], lam1[idx],   N=n.gwas,a=a, rough=rough)$a
    a=temp1
    h2 <- calH2.new1(x1, lam1,   N=n.gwas, a=a, rough=rough)$h2
  }else{
    result<- calH2.new1(x1, lam1,N=n.gwas, a=a, rough=rough)
    h2 <- as.numeric(result$h2)
    a <- as.numeric(result$a)
  }
  return(list(a=a,h2=h2))}

get.res.ldsc <- function(z1,ldsc1,n.gwas,a,twostage){
  m <- length(z1)
  s1thld <- 30
  if(twostage){
    a <-  calH2(z1[(z1^2)<s1thld], ldsc1[(z1^2)<s1thld],   N=n.gwas,a=NULL)$a
    h2 <- calH2(z1, ldsc1,   N=n.gwas,a=a)$h2
  }else{
    h2 <-  calH2(z1, ldsc1,   N=n.gwas,a=a)$h2
    a <- calH2(z1, ldsc1,   N=n.gwas,a=a)$a
  }
  return(list(a=a,h2=h2))}

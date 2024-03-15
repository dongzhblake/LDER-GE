ldscoreBand <- function(X, band=1){
  p <- ncol(X)
  N <- nrow(X)
  ldsc<- rep(0, p)
  for(i in 1:p){
    minBandInd<-max(1, i-band)
    maxBandInd<-min(p, i+band)
    if(minBandInd<i){
      Rl<-cor(X[,i], X[,minBandInd:(i-1)])
      Rl2 <- Rl^2
      Rl2 <- Rl2 - (1-Rl2)/(N-2)
    }else{
      Rl2 <-0
    }
    if(maxBandInd>i){
      Ru<-cor(X[,i], X[, (i+1):maxBandInd])
      Ru2 <- Ru^2
      Ru2 <- Ru2 - (1-Ru2)/(N-2)
    }else{
      Ru2 <-0
    }
    ldsc[i] <- 1+sum(Rl2)+sum(Ru2)
  }
  return(ldsc)
}

colSds <- function(X, na.rm=TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(X)) # thanks @flodel
  } else {
    n <- nrow(X)
  }
  colVar <- colMeans(X*X, na.rm=na.rm) - (colMeans(X, na.rm=na.rm))^2
  #	return(sqrt(colVar * n/(n-1)))
  return(sqrt(colVar))
}

normalize<-function(X){
  n <- nrow(X)
  p <- ncol(X)
  meanX <- rep(1,n)%o%colMeans(X)
  sdX <- rep(1, n)%o%colSds(X)
  normX <- (X-meanX)/sdX
  return(normX)
}

#Marginal Regression
margReg<-function(X, y){
  n <- length(y)
  normX <- normalize(X)
  normY <- normalize(y)
  beta <- t(normX)%*%normY/n
  se <- sqrt((1-beta^2)/(n-2))
  z <- beta/se
  return(z)
}
margReg.beta<-function(X, y){
  n <- length(y)
  normX <- normalize(X)
  normY <- normalize(y)
  beta <- t(normX)%*%normY/n
  return(sqrt(n)*beta)
}
ldscore<-function(R, N){
  R2 <- R^2
  R2 <- R2-(1-R2)/(N-2)
  #	diag(R2) <- 1
  if(!is.null(dim(R))){
    ldsc <- rowSums(R2)}else{
      ldsc <- R2
    }
  return(ldsc)
}

cov2corr<-function(A){
  sqrtDiagA<-sqrt(diag(A))
  return(A/(sqrtDiagA%o%sqrtDiagA))
}

weighting<-function(x, w){
  w <- sqrt(w)
  w <- w/sum(w)
  if(is.matrix(x)){
    if(nrow(x)!=length(w)){
      stop('ldsr.weighting: Row number in x does not match length of w!')
    }else{
      w<-matrix(rep(w, ncol(x)), ncol=ncol(x))
    }
  }else{
    if(length(x)!=length(w)){
      stop('ldsr.weighting: Lengths of x and w are not matched!')
    }
  }
  return(x*w)
}

wls<-function(x,y,w=rep(1, nrow(x))){
  x_new <- weighting(x, w)
  y_new <- weighting(y, w)
  lmObj <- lm.fit(x_new, y_new)
  return(lmObj$coefficients)
}

irwls<-function(x,y,updateFunc,w=rep(1, nrow(x)), iter=5){
  for(i in 1:iter){
    coeff <- wls(x,y,w) ## coeff[1]: coeff of ldsc; [2]:coeff of 1
    w <- updateFunc(coeff)
  }
  return(coeff)
}

calH2.new1 = function(z, ldsc, N, a=NULL, rough=F){
  if(length(z)!=length(ldsc)){
    stop('ldsr.calH2:lengths of z and ldsc are not matched!')}
  M <- length(z)
  chi2 <- z^2
  updateFunc <- function(coeff, w){
    temp1 <- min(max(0,coeff[1]),N/M)
    if(!is.null(a)){temp2 <- a*N}
    else{temp2 <- max(1,coeff[2])-1}
    w1 <- 1/(1+temp2+temp1*ldsc)^2
    w2 <- pmin(1.2,ldsc)*rough + pmin(1,ldsc)*(rough==F)
    return(w1*w2)}
  
  if(!is.null(a)){
    int <- N*a+1.
    y <- chi2-int
    coeff <- irwls(as.matrix(ldsc), y, updateFunc)
  }else{
    x <- cbind(ldsc, rep(1, length(ldsc)))
    y <- chi2
    coeff <- irwls(x, y, updateFunc)
    a <- (coeff[[2]]-1)/N}
  h2est <- coeff[1]*M/N
  return(list(h2=h2est, a=a))}


calH2<-function(z, ldsc, N, a=NULL){
  if(length(z)!=length(ldsc)){
    stop('ldsr.calH2:lengths of z and ldsc are not matched!')
  }
  M <- length(z)
  chi2 <- z^2
  updateFunc <- function(coeff, w){
    temp1 <- min(max(0,coeff[1]),N/M)
    if(!is.null(a)){temp2 <- a*N}
    else{temp2 <- max(1,coeff[2])-1}
    w1 <- 1/(1+temp2+temp1*ldsc)^2
    w2 <- 1/pmax(ldsc, 1)
    return(w1*w2)
  }
  if(!is.null(a)){
    int <- N*a+1.
    y <- chi2-int
    coeff <- irwls(as.matrix(ldsc), y, updateFunc)
  }else{
    x <- cbind(ldsc, rep(1, length(ldsc)))
    y <- chi2
    coeff <- irwls(x, y, updateFunc)
    a <- (coeff[[2]]-1)/N
  }

  h2est <- coeff[[1]]*M/N
  return(list(h2=h2est, a=a))
}

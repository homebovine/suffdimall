simu <- function(n, nsim, p, beta, std){
#    x <- rmvl(n , rep(0, p), diag(std^2, p, p))
    x <- matrix(rnorm(n * p, 0, std), n, p)
    xb <- x %*% beta
    y <- rnorm(n, sin(2  * xb) + 2 * exp( 2 + xb), sqrt(log( 2 + xb ^2)))
    mxy <- cbind(y, x)
    return(mxy)
    
    
}
simu1 <- function(n, p, beta, std){
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std ^2
    xb <-  rlogis(n, 0, b)
    invbb <- solve(t(beta) %*% beta)
    sigma <- (diag(1, p, p) - beta %*% invbb %*% t(beta))* std^2
    x <- matrix(NA, n, p)
    for(i in 1:n){
        x[i, ] <- mvrnorm (1, beta * as.numeric(invbb) * xb[i], sigma)
}
   #     u <- runif(p, -sqrt(3), sqrt(3) )
    #    bu <- sigma %*% u
    #x[i, ] <- beta * as.numeric(invbb) * xb[i] + bu

    y <- rnorm(n, sin(2  * xb) + 2 * exp( 2 + xb), sqrt(log( 2 + xb ^2)))
    mxy <- cbind(y, x)
    return(mxy)
}

simu2 <- function(n, p, beta, std, k){
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std^2
    betaest <- beta[2:p]
    xb <-  rlogis(n, 0, b)
    invbb <- solve(t(beta) %*% beta)
    sigma <- (diag(1, p, p) - beta %*% invbb %*% t(beta))* std^2
    lowersigma <- (diag(1, p -1, p-1) - betaest %*% invbb %*% t(betaest))* std^2
    if(p > 2){
    temp <- svd(lowersigma)
    D <- temp$u %*% (diag(temp$d)^(1/2)) %*% t(temp$v)
}else{
    D <- sqrt(sigma[2, 2])
}
    x <- matrix(NA, n, p)

    for(i in 1:n){
            m <- -digamma(1)
            v <- pi^2 /6 
         u <-   - log(-log(runif(p)))
        #u <- rgamma(p, shape = k, scale = 1)
       
        x[i, 2:p] <- betaest * as.numeric(invbb * xb[i]) +  D %*% (v^(-1/2) * (u[2:p] - m))
        x[i, 1] <- (xb[i] - betaest %*%  x[i, 2:p])/fix
}
   #     u <- runif(p, -sqrt(3), sqrt(3) )
    #    bu <- sigma %*% u
    #x[i, ] <- beta * as.numeric(invbb) * xb[i] + bu

    y <- rnorm(n, sin(2  * xb) + 2 * exp( 2 + xb), sqrt(log( 2 + xb ^2)))
    mxy <- cbind(y, x)
    return(mxy)
}

seff <- function(betaest, x, y){

    beta <- c(fix, betaest)
    p1 <- length(betaest)
    p <- p1 + 1
    n <- nrow(x)
    seff <- matrix(NA, n, p-1)
     for(i in 1: n){
        x2 = x[i, 2:p]
        xb = sum(x[i, ] * beta)
        invbb <- solve(t(beta) %*% beta)
        mx <- invbb * betaest * xb
        er = x2 - mx
        lf3 =  xb * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i]) ^2 /((xb ^2 + 2) * (log(xb^2 + 2))^2) - (-2 * cos(2 * xb) - 2 * exp(xb + 2)) * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i])/(log(xb^2 + 2))
        seff[i, ] <- er * lf3
        
    }
     mseff<- sum(apply(seff, 2, sum)^2)
    return(mseff)
}

seff3 <- function(betaest, x, y, std, norm){
    beta <- c(fix, betaest)
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std ^2
    rate = (t(beta) %*% beta)^(-1/2)
    p1 <- length(betaest)
    p <- p1 + 1
    n <- nrow(x)
    TT <- mer <- seff3 <- matrix(NA, n, p - 1)
    mf1 <- rep(NA, n)
    for(i in 1: n){
        x2 = x[i, 2:p]
        xb = sum(x[i, ] * beta)
        invbb <- solve(t(beta) %*% beta)
        mx <- invbb * betaest * xb
        er = x2 - mx
        if(norm == 0){
        cond = t(xb) %*% invbb %*% xb
        dcond = 2  * invbb %*% xb
        conv =  std ^2 * (diag(1, p-1, p-1) -  betaest %*%invbb %*% t(betaest))
        
        sigma12 = ing( p, 1, cond, conv, dcond)
        sigma = sigma12[[1]]
        dsigma = sigma12[[2]]
        invsigm = ginv(sigma)
    }else{
         csigm =   (diag(1, p - 1, p - 1) - betaest %*%invbb %*% t(betaest))
         sigm <- std ^2 * csigm
         invcsigm = ginv(csigm)
          invsigm <- std^(-2) * invcsigm
    }
        dmdxb = betaest * invbb
        lf1 = (1 - exp(xb/b))/(b * exp(xb/b) + b) # -std ^(-2) * invbb * xb
        lf3 =  xb * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i]) ^2 /((xb ^2 + 2) * (log(xb^2 + 2))^2) - (-2 * cos(2 * xb) - 2 * exp(xb + 2)) * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i])/(log(xb^2 + 2))
        dmdb2 = as.numeric(xb * invbb) * diag(1, p1, p1) - 2 * as.numeric(xb * invbb^2 )* betaest %*% t(betaest)# + x2 %*% invbb %*% t(betaest)
        if(norm == 0){
        seff3[i, ] <- er * lf1 + mx %*% t(er) %*% invsigm %*% dmdxb + er * lf3 + dmdb2 %*% invsigm %*% er +  dsigma   %*% invsigm %*% er
    }else
         seff3[i, ] <- er * lf1 + mx %*% t(er) %*% invsigm %*% dmdxb + er * lf3 + dmdb2 %*% invsigm %*% er
        mer[i, ]<- er
        mf1[i] <- lf1
        TT[i, ] <- mx %*% t(er) %*% invsigm %*% dmdxb  + dmdb2 %*% invsigm %*% er
    }
    mseff<- sum(apply(seff3, 2, sum)^2)
    return(mseff)
    
    
}

seff4 <- function(beta, x, y,std){
  
    p <- length(beta)
    n <- nrow(x)
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std ^2
    seff4 <- matrix(NA, n, p )
    tbb <- t(beta) %*% beta
    btb <- beta %*% t(beta)
     invbb <- solve(tbb)
     csigm =   (diag(1, p, p) - beta %*%invbb %*% t(beta))
    sigm <- std ^2 * csigm
    invcsigm = ginv(csigm)
    invsigm <- std^(-2) * invcsigm
    dmdxb = beta * invbb
     Dbb = std * csigm##p by p
   
    invDbb = std^(-1) * invcsigm ##p by p
     dmdxb = beta * invbb
       mC <- matrix(NA,p^2, p)
    for(i in 1: (p)){
        temp <- beta
        temp[i] <- 2 * beta[i]
        centerm <- matrix(0, p, p)
        centerm[i, ] <- temp
        centerm[, i]<- temp
 
        temp1 <- as.numeric(invbb) * centerm - btb * 2*  beta[i] * as.numeric(invbb^(2))
        uppD<- std * t(Dbb) %*% invDbb %*% temp1 %*% invDbb
        mC[, i]  <- as.vector((uppD))
        
    }
    computv<- function(xi, l){
        a <- xi
        b <- xi

        as.vector(a %*% t(b) - diag(1, l, l))
    }
    computevv<- function(i, a, b){

        (a[i, ] %*% t(b[i, ]))
    }
     computevvv<- function(i, a, b, c){

        as.vector((a[i, ] %*% t(b[i, ]))) %*% t(c[i, ])
    }
   # en = 1000
  #  simxmx <- (mvrnorm(en, rep(0, p), sigm)) ##n by p
 #   xi <-t(invDbb %*% t(simxmx)) ##n by p
    #dlf2 <-  t(invbb %*% t(beta) %*% invsigm %*% t(simxmx)) ## n by 1
#    dlf2e <- - xi
 #   simv<- t(apply(xi, 1, computv, p)) #n by p^2
    #Dlog2 = simxmx * matrix(dlf2, nrow = en, ncol = p) ##n by p
 #   Dlog2e = lapply(1 :en, computevvv, xi, dlf2e, simv)
  #  Dlog2v <- lapply(1 :en, computevvv, simxmx, dlf2, simv) ## n by p * (p)^2
  #  simvv <- lapply(1 :en, computevv, simv, simv)## n by p^2 * p^2
 #   mDlog2v <- Reduce("+", Dlog2v)/en
#    mDlog2e <-  Reduce("+", Dlog2e)/en
  #  invmsimvv <- ginv(Reduce("+", simvv)/en)
  #  k1 = mDlog2v %*% invmsimvv
  #  k1[, ] = 0
    k3 = 0
    k4 = t(mC) - t(dmdxb) %*% t(invDbb) %x% Dbb
    k5 = k4#-k4%*%  mDlog2e %*% invmsimvv
    for(i in 1: n){
    xb = sum(x[i, ] * beta)
   
    mx <- invbb * beta * xb
    er = x[i, ] - mx
    
    lf1 = (1 - exp(xb/b))/(b * exp(xb/b) + b)
    lf3 = xb * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i]) ^2 /((xb ^2 + 2) * (log(xb^2 + 2))^2) - (-2 * cos(2 * xb) - 2 * exp(xb + 2)) * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i])/(log(xb^2 + 2))
    
    normer = invDbb %*% er
    v <-as.vector( normer %*% t(normer) - diag(1, p, p))
    
    dmdb2 = as.numeric(xb * invbb) * diag(1, p, p) - 2 * as.numeric(xb * invbb^2 )* beta %*% t(beta) #+ x[i, ] %*% invbb %*% t(beta)
  
   
   k2 = - (t(dmdxb) %x%  mx + dmdb2)  %*% t(invDbb) #p-d by p
    
   
   
    seff4[i, ] <- er %*% lf1 + er * lf3 - k2 %*% normer  - k5 %*% v
}
    mseff<- sum(apply(seff4, 2, mean)^2)
    return(mseff)
}




seff41 <- function(betaest, x, y,std){

    p <- length(betaest) + 1
    beta <- c(fix, betaest)
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std ^2
    rate = (t(beta) %*% beta)^(-1/2)
    n <- nrow(x)
    seff4 <- matrix(NA, n, p-1 )
    tbb <- t(beta) %*% beta
    btb <- betaest %*% t(betaest)
     invbb <- solve(tbb)
     csigm =   (diag(1, p -1, p-1) - betaest %*%invbb %*% t(betaest))
    svdcsigm = svd(csigm)
  
    sigm <- std ^2 * csigm
      if(p >2){
    D = std * svdcsigm$u %*% diag(svdcsigm$d)^(1/2) %*% t(svdcsigm$v)
}else{
    D = sigm^(1/2)
}
    invcsigm = ginv(csigm)
    invsigm <- std^(-2) * invcsigm
    dmdxb = betaest * invbb
    invD <- ginv(D)
    # svdsigm <- svd(invsigm)
    #invD = svdsigm$u %*% diag(svdsigm$d)^(1/2) %*% t(svdsigm$v)
    # D <- solve(invD) #svdsigm$u %*% diag(svdsigm$d)^(-1/2) %*% t(svdsigm$v)
       mC <- matrix(NA,(p-1)^2, p -1 )
    for(i in 1: (p - 1)){
        temp <- temp0 <- betaest
        temp[i] <-  betaest[i] + eps
        betatemp <- c(fix, temp)
         tbbtemp <- t(betatemp) %*% betatemp
        invbbtemp <- ginv(tbbtemp)
          csigm1 =   (diag(1, p -1, p-1) - temp %*%invbbtemp %*% t(temp))
         svdcsigm1 = svd(csigm1)
        if(p > 2){
        D1 = std * svdcsigm1$u %*% diag(svdcsigm1$d)^(1/2) %*% t(svdcsigm1$v)
    }else{
        D1 = csigma1 * std
    }
    
        invD1 = ginv(D1)
        DD = (invD1 - invD)/eps
 
       
        uppD<- t(D) %*% DD
        mC[, i]  <- as.vector((uppD))
        
    }
   
    
    k4 = t(mC) - t(dmdxb) %*% t(invD) %x% D
    k5 = k4
    for(i in 1: n){
    xb = sum(x[i, ] * beta)
   
    mx <- invbb * betaest * xb
    er = x[i, -1] - mx
    
     lf1 = (1 - exp(xb/b))/(b * exp(xb/b) + b) # -rate ^2 * exp(-rate - 1)
    lf3 = xb * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i]) ^2 /((xb ^2 + 2) * (log(xb^2 + 2))^2) - (-2 * cos(2 * xb) - 2 * exp(xb + 2)) * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i])/(log(xb^2 + 2))
    
    normer = invD %*% er
    v <-as.vector( normer %*% t(normer) - diag(1, p - 1, p-1))
    
    dmdb2 = as.numeric(xb * invbb) * diag(1, p-1, p-1) - 2 * as.numeric(xb * invbb^2 )* betaest %*% t(betaest) #+ x[i, ] %*% invbb %*% t(beta)
  
   
   k2 = - (t(dmdxb) %x%  mx + dmdb2)  %*% t(invD) %*% normer #p-d by p
  #  k2 =-( mx %*% t(er)%*% invsigm %*% dmdxb + dmdb2 %*% t(invD) %*% normer)
   
  
    
    seff4[i, ] <- er %*% lf1 + er * lf3 - k2  - k5 %*% v
}
    mseff<- sum(apply(seff4, 2, sum)^2)
    return(mseff)
}

gthirdm <- function(theta, k){
    theta^3 * k * ( k+ 1 ) * (k+2)
    
}
gfourm <- function(theta, k){
    theta^4 * k * ( k+ 1 ) * (k+2) * (k + 3)
    
}
gsecm <- function(theta, k){
    theta^2 * k * ( k+ 1 ) 
    
}
fmm <- function(theta, k){
 v <- k * theta ^2
 u3 <- gthirdm(theta, k)
 u2 <- gsecm(theta, k)
 u <- k * theta
 u4 <- gfourm(theta, k)
 m <- u
 ru2 <- v^(-1) * (u2 - 2 * u * m + m ^2)
 ru4 <- v^(-2) * (u4 - 4 * m * u3 + 6 * m ^2 * u2 - 4 * m ^3 * u + m ^4)
 ru3 <- v^(-3/2) * (u3 - 3 * u2 * m + 3 * u * m ^2 - m ^3)
 return(c(ru2, ru3, ru4))
    
}
seff42 <- function(betaest, x, y,std, k){

    p <- length(betaest) + 1
    beta <- c(fix, betaest)
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std ^2
    n <- nrow(x)
    seff4 <- matrix(NA, n, p-1 )
    tbb <- t(beta) %*% beta
    btb <- betaest %*% t(betaest)
     invbb <- solve(tbb)
     csigm =   (diag(1, p -1, p-1) - betaest %*%invbb %*% t(betaest))
    svdcsigm = svd(csigm)
    
 
    sigm <- std ^2 * csigm
      if(p >2){
    D = std * svdcsigm$u %*% diag(svdcsigm$d)^(1/2) %*% t(svdcsigm$v)
}else{
    D = sigm^(1/2)
}
    invcsigm = ginv(csigm)
    invsigm <- std^(-2) * invcsigm
    dmdxb = betaest * invbb
    invD <- ginv(D)
    # svdsigm <- svd(invsigm)
    #invD = svdsigm$u %*% diag(svdsigm$d)^(1/2) %*% t(svdsigm$v)
    # D <- solve(invD) #svdsigm$u %*% diag(svdsigm$d)^(-1/2) %*% t(svdsigm$v)
       mC <- matrix(NA,(p-1)^2, p -1 )
   xb <- 1
     mm3 <- fmm(1, k)
    secm <- 1# mm3[1]
    thirdm <- 1.14#mm3[2]
    fourm <- 12/5 + 3#mm3[3]
    v1em <- vector('list')
    v1v1m <- vector('list')
    eem <- diag(secm, p-1, p-1)
    iden <- diag(1, p-1, p-1)
    for(j in 1: (p -1)){
        v1em[[j]] <- matrix(0, (p - 1), (p - 1))
        v1em[[j]][j, j] <- thirdm
        v1v1m[[j]] <- matrix(0, (p - 1), (p - 1))
        v1v1m[[j]][j, j] <- fourm
        if( p > 2){
        diag(v1v1m[[j]]) <- (secm) ^2
         v1v1m[[j]][j, j] <- fourm
    }
        
    }
    v1em <- do.call(rbind, v1em) #(p-1)^2 by p-1
    v1v1m <- as.matrix(bdiag(v1v1m))
     Evv <- v1v1m - as.vector(iden) %*% t(as.vector(iden))-   v1em %*% t(v1em) 
    invEvv <- solve(Evv)
    for(i in 1: n){
    xb = sum(x[i, ] * beta)
  # (p-1)^2 by (p-1)^2
   
    mx <- invbb * betaest * xb
    er = x[i, -1] - mx
    
     lf1 = (1 - exp(xb/b))/(b * exp(xb/b) + b) # -rate ^2 * exp(-rate - 1)
    lf3 = xb * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i]) ^2 /((xb ^2 + 2) * (log(xb^2 + 2))^2) - (-2 * cos(2 * xb) - 2 * exp(xb + 2)) * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i])/(log(xb^2 + 2))
    
    normer = invD %*% er
    v <- as.vector( normer %*% t(normer)  - diag(1, p - 1, p-1)) - v1em %*% normer
    
    dmdb2 = as.numeric(xb * invbb) * diag(1, p-1, p-1) - 2 * as.numeric(xb * invbb^2 )* betaest %*% t(betaest) #+ x[i, ] %*% invbb %*% t(beta)
  
   
   k2 = - (t(dmdxb) %x%  mx + dmdb2)  %*% t(invD) #p-d by p
  #  k2 =-( mx %*% t(er)%*% invsigm %*% dmdxb + dmdb2 %*% t(invD) %*% normer)
   
  k3 <- k2 %*% t(v1em) %*% as.matrix(invEvv)
    
    seff4[i, ] <- er %*% lf1 + er * lf3 - k2   %*% normer + k3 %*% v
}
    mseff<- sum(apply(seff4, 2, sum)^2)
    return(mseff)
}



ing <- function( p, d, cond, conv, dcond){
    f1 <- function(t, cond, p, d) exp(-(t  + cond)^(1/2) ) * t ^((p-d)/2)
    f2 <- function(t, cond, p, d) exp( -(t  + cond)^(1/2) ) * t ^((p-d)/2 -1)
    f3 <- function(t, cond, p, d, dcond) exp( - (t  + cond)^(1/2) ) * t ^((p-d)/2)  * (-1/2) * (t + cond)^(-1/2) * dcond
    f4 <- function(t, cond, p, d, dcond) exp( -(t  + cond)^(1/2) ) * t ^((p-d)/2 - 1)  * (-1/4) * (t + cond)^(-1/2) * dcond
    num<- integrate(f1, 0, Inf, cond, p, d)$value
    denom <- integrate(f2, 0, Inf, cond, p, d)$value
    dnum<- integrate(f3, 0, Inf, cond, p, d, dcond)$value
    ddenom <- integrate(f4, 0, Inf, cond, p, d, dcond)$value
   
    res <- num/denom * conv * 1/(p -d)
    dres <- (dnum*  denom - ddenom * num)/denom^2 *  conv * 1/(p -d)
    list(res, dres)
}
monting <- function( p, d, cond, conv, dcond){
    f1 <- function(t, cond, p, d) exp(-1/2 * (t  + cond)^(1/2) ) * t ^((p-d)/2)
    f2 <- function(t, cond, p, d) exp(-1/2 * (t  + cond)^(1/2) ) * t ^((p-d)/2 -1)
    f3 <- function(t, cond, p, d, dcond) exp(-1/2 * (t  + cond)^(1/2) ) * t ^((p-d)/2)  * (-1/4) * (t + cond)^(-1/2) * dcond
    f4 <- function(t, cond, p, d, dcond) exp(-1/2 * (t  + cond)^(1/2) ) * t ^((p-d)/2 - 1)  * (-1/4) * (t + cond)^(-1/2) * dcond
    vt<- runif(100, 0, 1e4)
    num <- mean(f1(vt, cond, p, d)) * 1e4
    
    denom <- mean(f2(vt, cond, p, d)) * 1e4
    dnum<-mean(f3(vt, cond, p, d, dcond)) * 1e4 
    ddenom <-mean(f4(vt, cond, p, d, dcond)) * 1e4
   
    res <- num/denom * conv * 1/(p -d)
    dres <- (dnum*  denom - ddenom * num)/denom^2 *  conv * 1/(p -d)
    list(res, dres)
}
rmvl <- function(n, mu, Sigma)
     {
     mu <- rbind(mu)
     if(missing(Sigma)) Sigma <- diag(ncol(mu))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     
     k <- ncol(Sigma)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     e <- matrix(rexp(n, 1), n, k)
     z <- rmvn(n, rep(0, k), Sigma)
     x <- mu + sqrt(e)*z
     return(x)
     }

rmvn <- function(n=1, mu=rep(0,k), Sigma)
     {
     mu <- rbind(mu)
     if(missing(Sigma)) Sigma <- diag(ncol(mu))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
   
     k <- ncol(Sigma)
     if(n > nrow(mu)) mu <- matrix(mu, n, k, byrow=TRUE)
     z <- matrix(rnorm(n*k),n,k) %*% chol(Sigma)
     x <- mu + z
     return(x)
     o}

estimation <- function(beta0, betaini, n, p, itr){
    set.seed(itr + 3016)
    
    simdata <- simu2(n,  p, beta0, std, 1/2)
    y <- simdata[, 1]
    x <- simdata[, 2: (p+ 1)]
    temp1 <-optim(betaini, seff, gr = NULL, x, y, method = 'Nelder-Mead')# try(nlm( seff,betaini,  x, y))### #spg(betaini, seff, gr=NULL, method=3, lower=-Inf, upper=Inf,  project=NULL, projectArgs=NULL,  control=list(), quiet=FALSE, alertConvergence=TRUE, x, y)$par##dfsane(betaini, seff, method=2, control=list(trace = FALSE, tol= 1e-4),quiet=FALSE, alertConvergence=TRUE, x, y)$par###
    #estseff <- temp1$estimate
   # temp1 <-try(nlm( seff,estseff,  x, y))# optim(c( estseff), seff, gr = NULL, x, y, method = 'BFGS')#
    estseff <- temp1$par
    temp1 <-optim(c( estseff)* 1.01, seff, gr = NULL, x, y, method = 'BFGS')
    estseff1 <- temp1$par
    temp2 <-optim(c( betaini), seff3, gr = NULL, x, y, std, 1, method = 'Nelder-Mead')
    temp2 <-optim(temp2$par * 1.01, seff3, gr = NULL, x, y, std, 1, method = 'BFGS')
                                        #nlm(seff3,betaini, x, y, std, 1)# ######spg(betaini, seff3, gr=NULL, method=3, lower=-Inf, upper=Inf,  project=NULL, projectArgs=NULL,  control=list(), quiet=FALSE, alertConvergence=TRUE, x, y, std)$par#dfsane(betaini, seff3, method=2, control=list(trace = FALSE, tol = 1e-4), quiet=FALSE, alertConvergence=TRUE, x, y, std)$par#
    estseff3 <- temp2$par
    temp3 <- optim(c( betaini), seff42, gr = NULL, x, y, std, 1/2,  method = 'Nelder-Mead')
    temp3 <- optim(temp3$par * 1.01, seff42, gr = NULL, x, y, std, 1/2, method = 'BFGS')#nlm(seff41,c(betaini), x, y, std)####optim(c(1, estseff), seff4, gr = NULL, x, y, std, method = 'L-BFGS-B')##  spg(c(1, betaini), seff4, gr=NULL, method=3, lower=-Inf, upper=Inf,  project=NULL, projectArgs=NULL,  control=list(), quiet=FALSE, alertConvergence=TRUE, x, y, std)$par#, method = 'L-BFGS-B', lower = c(1, rep(-Inf, p-1)), upper = c(1, rep(Inf, p-1)))$par
    estseff4 <- temp3$par
   # estseff4 <- estseff4[-1]/estseff4[1]
    print(c(temp1$convergence, temp2$convergence,  temp3$convergence))
    print(c(temp1$value, temp2$value, temp3$value))#, temp3$convergence))
    #estseff4 <- (estseff4/estseff4[1])[-1]
    return(c(estseff1,estseff3,  estseff4))#,  estseff4))
}
nsim <- 1200
p <- 4
res <- matrix(NA, nsim, (p - 1) * 2)
res1 <- matrix(NA, nsim, (p-1) * 2)
res2 <- matrix(NA, nsim, (p-1) * 2)
res3 <- matrix(NA, nsim, (p-1) * 3)
res3 <- rbind(res3, res1)
fix <-0.3
for(itr in 1:nsim){
    print(itr)
    std <- 0.1
    n <-50
    
    beta0 <- c(fix, 0.2, 0.3, 0.4)
    betaini <- c( beta0[2:p] * 1.21)
    
    temp<- try(estimation(beta0, betaini, n, p, itr))
    if(class(temp) != 'try-error')
        res3[itr, ] <- temp
}
ix <-  500
mbeta0 <- matrix(beta0[2:4], ix, ncol = p-1, byrow = T)
sqereff <- sqrt(apply((res3[1:ix, 1:3] - mbeta0)^2, 1, mean, na.rm = T))
sqereff3 <- sqrt(apply((res3[1:ix, 4:6] - mbeta0)^2, 1, mean, na.rm = T))
sqereff4 <- sqrt(apply((res3[1:ix, 7:9] - mbeta0)^2, 1, mean, na.rm = T))#apply((res[, 11:15] - mbeta0)^2, 1, sum)
sqereff <- sqereff[which(sqereff<= quantile(sqereff, 0.95, na.rm = T))]
sqereff3 <- sqereff3[which(sqereff3<= quantile(sqereff3, 0.95, na.rm = T))]
sqereff4 <- sqereff4[which(sqereff4<= quantile(sqereff4, 0.95, na.rm = T))]
pdf('temp.pdf')
boxplot(sqereff, sqereff3, sqereff4, range = 1, outline = F, names = c("a", "b", "c"))
dev.off()

library(MASS)
library(Matrix)

simu <- function(n, nsim, p, beta, std){
    x <- rmvl(n , rep(0, p), diag(std^2, p, p))
#    x <- matrix(rnorm(n * p, 0, std), n, p)
    xb <- x %*% beta
    y <- rnorm(n, sin(2  * xb) + 2 * exp( 2 + xb), sqrt(log( 2 + xb ^2)))
    mxy <- cbind(y, x)
    return(mxy)
    
    
}
simu1 <- function(n, p, beta, std){
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std
    xb <-  rlogis(n, 0, b)
    invbb <- solve(t(beta) %*% beta)
    sigma <-  (diag(1, p, p) - beta %*% invbb %*% t(beta))* std^2
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

simu2 <- function(n, p, beta, std, std1,  k, nonconst, linear, t0, t1){
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std1^2
    betaest <- beta[2:p]
    xb <-  rlogis(n, 0, b)
      k0 = t(beta) %*% beta/t0^2
    kt = k0 * t0
  #  xb = rgamma(n, k0, scale = t0) - kt
    invbb <- solve(t(beta) %*% beta)
    
    sigma <- (diag(1, p, p) - beta %*% invbb %*% t(beta))* std^2
    if(nonconst == 0){
    lowersigma <- (diag(1, p -1, p-1) - betaest %*% invbb %*% t(betaest))* std^2
    if(p > 2){
    temp <- svd(lowersigma)
    D <- temp$u %*% (diag(temp$d)^(1/2)) %*% t(temp$v)
}else{
    D <- sqrt(sigma[2, 2])
}
}
    x <- matrix(NA, n, p)

    for(i in 1:n){
        if(nonconst == 1){
            lowersigma <- as.numeric(xb[i] ^2/(t(beta) %*% beta)) * (diag(1, p -1, p-1) - betaest %*% invbb %*% t(betaest))
            if(p > 2){
             temp <- svd(lowersigma)
             D <- temp$u %*% (diag(temp$d)^(1/2)) %*% t(temp$v)
         }else{
             D <- sqrt(sigma[2, 2] *  xb[i] ^2/(t(beta) %*% beta))
         }
    }
            m <-t1  * gamma(1 + 1/k)#k  * t1# -digamma(1)
            v <-t1^2* (gamma(1 + 2/k) - (gamma(1 + 1/k))^2)#k  * t1^2# pi^2 /6 
         #u <-   - log(-log(runif(p)))
        u <- rweibull(p, shape = k, scale = t1)#rgamma(p, shape = k, scale = t1)
if(linear == 0){       
        x[i, 2:p] <- (betaest * as.numeric((invbb * xb[i]^2))) +  as.matrix(D) %*% as.matrix(v^(-1/2) * (u[2:p] - m))
        x[i, 1] <- (xb[i] - betaest %*%  x[i, 2:p])/fix
    }else{
           x[i, 2:p] <- betaest * as.numeric(invbb * xb[i]) +  as.matrix(D) %*% as.matrix(v^(-1/2) * (u[2:p] - m))
        x[i, 1] <- (xb[i] - betaest %*%  x[i, 2:p])/fix

    }
}
   #     u <- runif(p, -sqrt(3), sqrt(3) )
    #    bu <- sigma %*% u
    #x[i, ] <- beta * as.numeric(invbb) * xb[i] + bu

   y <- rnorm(n, sin(2  * xb) * c, d * sqrt(log( 2 + xb ^2)))
    mxy <- cbind(y, x)
    return(mxy)
}

simu3 <- function(n, p, beta, std, std1,  k,  t0, t1){
    k0 = t(beta) %*% beta/t0^2
    kt = k0 * t0
    xb = rgamma(n, k0, scale =  t0) - kt
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std1^2
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
        

     
    
            m <-t1  * gamma(1 + 1/k)#k  * t1# -digamma(1)
            v <-t1^2* (gamma(1 + 2/k) - (gamma(1 + 1/k))^2)#k  * t1^2# pi^2 /6 
         #u <-   - log(-log(runif(p)))
        u <- rweibull(p, shape = k, scale = t1)

        e2 = D %*% (v^(-1/2) * (u[2:p] - m))
        e1 = -t(betaest) %*% e2
           x[i, ] <- beta %*% (invbb) %*% xb[i] + c(e1, e2)
        

    
}
   #     u <- runif(p, -sqrt(3), sqrt(3) )
    #    bu <- sigma %*% u
    #x[i, ] <- beta * as.numeric(invbb) * xb[i] + bu

    y <- rnorm(n, sin(2  * xb) + 2 * exp( 2 + xb), sqrt(log( 2 + xb ^2)))
    mxy <- cbind(y, x)
    return(mxy)
}


seff <- function(betaest, x, y, linear){

    beta <- c(fix, betaest)
    p1 <- length(betaest)
    p <- p1 + 1
    n <- nrow(x)
    seff <- matrix(NA, n, p-1)
     for(i in 1: n){
        x2 = x[i, 2:p]
        xb = sum(x[i, ] * beta)
        invbb <- solve(t(beta) %*% beta)
        if(linear == 1){
        mx <- invbb * betaest * xb
    }else{
           mx <- invbb * betaest * xb^2
    }
        er = x2 - mx
          lf3 =  xb * (-sin(2 * xb)* c + y[i]) ^2 /(d * (xb ^2 + 2) * (log(xb^2 + 2))^2) + (2 * c * cos(2 * xb)) * (-sin(2 * xb) * c  + y[i])/(d * log(xb^2 + 2)) - xb /((xb^2+ 2) * log(xb^2 + 2))#xb * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i]) ^2 /((xb ^2 + 2) * (log(xb^2 + 2))^2) - (-2 * cos(2 * xb) - 2 * exp(xb + 2)) * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i])/(log(xb^2 + 2)) - xb /((xb^2+ 2) * log(xb^2 + 2))
        seff[i, ] <- er * lf3
        
    }
     mseff<- sum(apply(seff, 2, sum)^2)
    return(mseff)
}

seff3 <- function(betaest, x, y, std, std1,  nonconst, t0, t1){
    
    
    beta <- c(fix, betaest)
    k0 = t(beta) %*% beta/t0^2
    kt = k0 * t0
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std1 ^2
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
        if(nonconst == 1){
        csigm =  as.numeric(xb ^2/ (t(beta) %*% beta)) *  (diag(1, p - 1, p - 1) - betaest %*%invbb %*% t(betaest))
         sigm <- std ^2 * csigm
         invcsigm = ginv(csigm)
          invsigm <- std^(-2) * invcsigm
        dsigma <- as.numeric(2 * xb /(t(beta) %*% beta))* (diag(1, p - 1, p - 1) - betaest %*%invbb %*% t(betaest))
    }else{
         csigm =   (diag(1, p - 1, p - 1) - betaest %*%invbb %*% t(betaest))
         sigm <- std ^2 * csigm
         invcsigm = ginv(csigm)
          invsigm <- std^(-2) * invcsigm
    }
        dmdxb = betaest * invbb
        #lf1 = (xb > -kt) *((k0-1)/(xb + kt) - 1/t0)
         lf1 = (1 - exp(xb/b))/(b * exp(xb/b) + b) # -std ^(-2) * invbb * xb
   lf3 =  xb * (-sin(2 * xb)* c + y[i]) ^2 /(d * (xb ^2 + 2) * (log(xb^2 + 2))^2) + (2 * c * cos(2 * xb)) * (-sin(2 * xb) * c  + y[i])/(d * log(xb^2 + 2)) - xb /((xb^2+ 2) * log(xb^2 + 2))
       dmdb2 = as.numeric(xb * invbb) * diag(1, p1, p1) - 2 * as.numeric(xb * invbb^2 )* betaest %*% t(betaest)# + x2 %*% invbb %*% t(betaest) (1 - betaest^2)/(1 + betaest^2)^2 * xb
        if(nonconst == 1){
        seff3[i, ] <- er * lf1 + mx %*% t(er) %*% invsigm %*% dmdxb + er * lf3 + dmdb2 %*% invsigm %*% er +  dsigma   %*% invsigm %*% er
    }else{
         seff3[i, ] <- er * lf1 + mx %*% t(er) %*% invsigm %*% dmdxb + er * lf3 + dmdb2 %*% invsigm %*% er
     }
       
    }
    mseff<- sum(apply(seff3, 2, sum)^2)
    return(mseff)
    
    
}


seff42 <- function(betaest, x, y,std, std1, k, t0, t1){
    
    m <-t1  * gamma(1 + 1/k)#k  * t1# -digamma(1)
    v <-t1^2* (gamma(1 + 2/k) - (gamma(1 + 1/k))^2)#k  * t1^2# pi^2 /6 
                
    p <- length(betaest) + 1
    beta <- c(fix, betaest)
    k0 = t(beta) %*% beta/t0^2
    kt = k0 * t0
    b <- sqrt(3 * t(beta) %*% beta /pi^2) * std1 ^2
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
        # mm3 <- fmm(1, k)
    secm <- 1# mm3[1]#1#
    thirdm <- (gamma(1 + 3/k) * t1^3 - 3 * m * v - m^3)/v^(3/2)#mm3[2]#1.14#mm3[2]
    gamma1 <- gamma(1+ 1/k)
     gamma2 <- gamma(1+ 2/k)
     gamma3 <- gamma(1+ 3/k)
     gamma4 <- gamma(1+ 4/k)
    fourm <- (-6 * gamma1^4 + 12 * gamma1^2 * gamma2- 3 * gamma2^2 - 4 * gamma1 * gamma3 + gamma4)/(gamma2 - gamma1^2)^2 + 3#mm3[3]#12/5 + 3#mm3[3]
    iden <- diag(1, p-1, p-1)
    if(p > 2){
   xb <- 1
    
   
    eem <- diag(secm, p-1, p-1)
    
    efv2 <- efv1 <- matrix(0, p-1, p-1)
    efv1[1, 2] <- -1
    efv1 [2, 1] <- -1
    efv2[1, 3] <- -1
    efv2[3, 1] <- -1
    
    nvI<- as.vector(diag(-1, p-1, p-1))
    efI <- matrix(0, (p - 1)^2, (p-1)^2)
    efI[, c(1, 5, 9)]<- nvI
    efe <- matrix(0, (p-1)^2, p-1)
    efe[1, 1] <- efe[5, 2] <- efe[9, 3] <- -1
}
    v1em <- vector('list')
    efv <- v1v1m <- vector('list')
    for(j in 1: (p -1)){
        if(p > 2){
        efv3 <- diag(-1, 3, 3)
        if(j == 1){
            efv3[j, j] <- -3
            efv[[j]] <- cbind(efv3, efv1, efv2)
        }else if(j == 2){
             efv3[j, j] <- -3
             efv[[j]] <- cbind( efv1, efv3,  efv2)
        }else{
             efv3[j, j] <- -3
             efv[[j]] <- cbind( efv1, efv2,  efv3)
        }
    }
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
    
    invEvv <- ginv(Evv)
    if(p > 2){
    efv <- do.call(rbind, efv)
}

if(p > 2){
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
        D1 = csigm1 * std
    }
    
        invD1 = ginv(D1)
        DD = (invD1 - invD)/eps
 
       
        uppD<- t(D) %*% DD
        mC[, i]  <- as.vector((uppD))
        
    }

    
    k4 = t(mC) - t(dmdxb) %*% t(invD) %x% D
    k5 <- -k4 %*% (efv - efI - efe %*% t(v1em)) %*% invEvv
}else k4 <- k5 <- 0 
    for(i in 1: n){
    xb = sum(x[i, ] * beta)
  # (p-1)^2 by (p-1)^2
   
    mx <- invbb * betaest * xb
    er = x[i, -1] - mx
   # lf1 = (xb > -kt) *((k0-1)/(xb + kt) - 1/t0)
  #  lf1 = - (xb/b) - 2 * log(1 + exp(-xb/b))
    lf1 = (1 - exp(xb/b))/(b * exp(xb/b) + b) # -rate ^2 * exp(-rate - 1)
          lf3 =  xb * (-sin(2 * xb)* c + y[i]) ^2 /(d * (xb ^2 + 2) * (log(xb^2 + 2))^2) + (2 * c * cos(2 * xb)) * (-sin(2 * xb) * c  + y[i])/(d * log(xb^2 + 2)) - xb /((xb^2+ 2) * log(xb^2 + 2))
    
    normer = invD %*% er
    v <- as.vector( normer %*% t(normer)  - diag(1, p - 1, p-1)) - v1em %*% normer
    
    dmdb2 = as.numeric(xb * invbb) * diag(1, p-1, p-1) - 2 * as.numeric(xb * invbb^2 )* betaest %*% t(betaest) #+ x[i, ] %*% invbb %*% t(beta)
  
   
   k2 = - (t(dmdxb) %x%  mx + dmdb2)  %*% t(invD) #p-d by p
  #  k2 =-( mx %*% t(er)%*% invsigm %*% dmdxb + dmdb2 %*% t(invD) %*% normer)
   
  k3 <- k2 %*% t(v1em) %*% as.matrix(invEvv)
    
    seff4[i, ] <- er %*% lf1 + er * lf3 - k2   %*% normer + k3 %*% v - k5 %*% v
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
    
    simdata <- simu2(n,  p, beta0, std, std1, k, 0, 1,   t0, t1)
    y <- simdata[, 1]
    x <- simdata[, 2: (p+ 1)]
    temp1 <-optim(betaini, seff, gr = NULL, x, y,1, method = 'Nelder-Mead')# try(nlm( seff,betaini,  x, y))### #spg(betaini, seff, gr=NULL, method=3, lower=-Inf, upper=Inf,  project=NULL, projectArgs=NULL,  control=list(), quiet=FALSE, alertConvergence=TRUE, x, y)$par##dfsane(betaini, seff, method=2, control=list(trace = FALSE, tol= 1e-4),quiet=FALSE, alertConvergence=TRUE, x, y)$par###
    estseff <- temp1$par
   # temp1 <-try(nlm( seff,estseff,  x, y))# optim(c( estseff), seff, gr = NULL, x, y, method = 'BFGS')#
  #  estseff <- temp1$par
  #  temp1 <-optim(betaini, seff, gr = NULL, x, 0, 1, method = 'BFGS')
 #   estseff1 <- temp1$par
    temp2 <-optim(c( betaini), seff3, gr = NULL, x, y, std, std1, 0, t0, t1,  method = 'Nelder-Mead')
#    temp2 <-optim(betaini, seff3, gr = NULL, x, y, std, std1, 0, t0, t1,   method = 'BFGS')
                                        #nlm(seff3,betaini, x, y, std, 1)# ######spg(betaini, seff3, gr=NULL, method=3, lower=-Inf, upper=Inf,  project=NULL, projectArgs=NULL,  control=list(), quiet=FALSE, alertConvergence=TRUE, x, y, std)$par#dfsane(betaini, seff3, method=2, control=list(trace = FALSE, tol = 1e-4), quiet=FALSE, alertConvergence=TRUE, x, y, std)$par#
    estseff3 <- temp2$par
    temp3 <- optim(c( betaini), seff42, gr = NULL, x, y, std, std1, k, t0, t1,  method = 'Nelder-Mead')
  #  temp3 <- optim(temp3$par * 1.01, seff42, gr = NULL, x, y, std, std1, k, t0, t1,   method = 'BFGS')#nlm(seff41,c(betaini), x, y, std)####optim(c(1, estseff), seff4, gr = NULL, x, y, std, method = 'L-BFGS-B')##  spg(c(1, betaini), seff4, gr=NULL, method=3, lower=-Inf, upper=Inf,  project=NULL, projectArgs=NULL,  control=list(), quiet=FALSE, alertConvergence=TRUE, x, y, std)$par#, method = 'L-BFGS-B', lower = c(1, rep(-Inf, p-1)), upper = c(1, rep(Inf, p-1)))$par
    estseff4 <- temp3$par
   # estseff4 <- estseff4[-1]/estseff4[1]
    print(c(temp3$convergence))
    print(c(temp1$par, temp2$par, temp3$par))#, temp3$convergence))
    #estseff4 <- (estseff4/estseff4[1])[-1]
    return(c(estseff, estseff3, estseff4))#,  estseff4))
}



nsim <- 1200
p <- 4
res <- matrix(NA, nsim, (p - 1) * 3)
res1 <- matrix(NA, nsim, (p-1) * 3)
res2 <- matrix(NA, nsim, (p-1) * 3)
t0 <- 5
t1 <- 2
#res <- rbind(res, res1)
fix <- 1
k <- 1
 for(itr in 295:nsim){
    print(itr)
    std <- 1
    std1 <- sqrt(1)
    n <- 50
    
    beta0 <- c(fix, 0.2, 0.3, 0.4)
    betaini <- c( beta0[2:p] * 1.21)
    
    temp<- try(estimation(beta0, betaini, n, p, itr))
    if(class(temp) != 'try-error')
        res[itr, ] <- temp
}
 ix <- 295
mbeta0 <- matrix(beta0[2:4], ix, ncol = p-1, byrow = T)
sqereff <- sqrt(apply((res[1:ix, 1:3] - mbeta0)^2, 1, mean, na.rm = T))
sqereff3 <- sqrt(apply((res[1:ix, 4:6] - mbeta0)^2, 1, mean, na.rm = T))
sqereff4 <- sqrt(apply((res[1:ix, 7:9] - mbeta0)^2, 1, mean, na.rm = T))#apply((res[, 11:15] - mbeta0)^2, 1, sum)
sqereff <- sqereff[which(sqereff<= quantile(sqereff, 0.95, na.rm = T))]
sqereff3 <- sqereff3[which(sqereff3<= quantile(sqereff3, 0.95, na.rm = T))]
sqereff4 <- sqereff4[which(sqereff4<= quantile(sqereff4, 0.95, na.rm = T))]
pdf('temp.pdf')
boxplot(as.vector(sqereff), as.vector(sqereff3), sqereff4,   range = 0.8, outline = F, names = c("a", "b", 'c'), pars = (list(boxwex = 0.5)))
dev.off()

resab <- res
res3 = cbind(resab, res)
res1 <- res

abb1 = round(abs(apply(res3[1:ix, 1:3], 2, median, na.rm = T) - beta0[2:p]), 4)
abb2 = round(abs(apply(res3[1:ix, 4:6], 2, median, na.rm = T) - beta0[2:p]), 4)
abb3 = round(abs(apply(res3[1:ix, 7:9], 2, median, na.rm = T) - beta0[2:p]), 4)

mad1 = round(apply(res3[1:ix, 1:3], 2, mad), 4)
mad2 = round(apply(res3[1:ix, 4:6], 2, mad), 4)
mad3 = round(apply(res3[1:ix, 7:9], 2, mad), 4)

print(noquote(paste(abb1, '(', mad1, ')','&',  abb2, '(', mad2, ')','&', abb3, '(', mad3, ')', '\\')))

upper.trim.mean <- function(x,trim) {
 
  mean(x[x< quantile(x, trim, na.rm = T)],  na.rm = T)
}

mvar <- function(fres){
  
    
    bias<- round(matrix(abs(apply(fres[1:1200, ], 2, mean, 0.1, na.rm = T) - beta0[2:p]), nrow = 3) , 5)
std <- round(matrix(apply(fres[], 2, mad, na.rm = T), nrow = 3), 4)
    list(bias, std)
   
}
mfres <- list(fres1, fres2, fres3, fres1gamma, fres2gamma, fres3gamma)
mbiasstd <- list()
for(i in 1:6 ){
    mbiasstd[[i]] <- mvar(mfres[[i]])
}
print(noquote(paste(abb1, '(', mad1, ')','&',  abb2, '(', mad2, ')','&', abb3, '(', mad3, ')', '\\')))
i <- 1
print(noquote(paste(mbiasstd[[i]][[1]][, 1], ' (', mbiasstd[[i]][[2]][, 1], ')', '&',mbiasstd[[i]][[1]][, 2], ' (', mbiasstd[[i]][[2]][, 2], ')', '&', mbiasstd[[i]][[1]][, 3], ' (', mbiasstd[[i]][[2]][, 3], ')', '\\',  sep = '' )))

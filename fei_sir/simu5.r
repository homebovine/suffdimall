library(MASS)
library(Matrix)



simu2 <- function(n, p, beta,  k, nonconst, linear, t0, t1){
  #  b <- sqrt(3 * t(beta) %*% beta /pi^2) * std1^2
    betaest <- beta[3:p, ] ## p-d by d
    #xb1 <-  rlogis(n, 0, b)
    a11 <- t(beta[, 1]) %*% beta[, 1]
    a12 <- t(beta[, 1]) %*% beta[, 2]
    a22 <- t(beta[, 2]) %*% beta[, 2]
    tbb <- t(beta) %*% beta
    tbbroot <- svd(tbb)
    tbbroot <- tbbroot$u %*% diag(tbbroot$d)^{1/2} %*% t(tbbroot$v)
    k0 = a11/t0^2
    kt = k0 * t0
    xb1 = rgamma(n, k0, scale = t0) - kt
    b <- sqrt(3 * (a22 - a12 * a11^{-1} * a12) /pi^2) 
    xb2 = rlogis(n, a11^{-1} * a12 * xb1, b)
    invbb <- solve(t(beta) %*% beta) #d by d
    xb = cbind(xb1, xb2)
    sigma <- (diag(1, p, p) - beta %*% invbb %*% t(beta))
    if(nonconst == 0){
    lowersigma <- (diag(1, p -2, p-2) - betaest %*% invbb %*% t(betaest)) #p-d by p-d
   
    temp <- svd(lowersigma)
    D <- temp$u %*% (diag(temp$d)^(1/2)) %*% t(temp$v)
}
    x <- matrix(NA, n, p)

    for(i in 1:n){
        if(nonconst == 1){
            lowersigma <- as.numeric(xb1[i] ^2/(t(beta[, 1]) %*% beta[, 1])) * (diag(1, p -2, p-2) - betaest %*% invbb %*% t(betaest))
           
             temp <- svd(lowersigma)
             D <- temp$u %*% (diag(temp$d)^(1/2)) %*% t(temp$v)
        
    }
            m <-t1  * gamma(1 + 1/k)#k  * t1# -digamma(1)
            v <-t1^2* (gamma(1 + 2/k) - (gamma(1 + 1/k))^2)#k  * t1^2# pi^2 /6 
         #u <-   - log(-log(runif(p)))
        u <- rweibull(p, shape = k, scale = t1)#rgamma(p, shape = k, scale = t1)
if(linear == 0){       
        x[i, 3:p] <- (betaest %*% invbb  %*% (xb[i, ]^2 - diag(tbb))) +  as.matrix(D) %*% as.matrix(v^(-1/2) * (u[3:p] - m))
        x[i, 1] <- (xb[i, 1] - betaest[, 1] %*%  x[i, 3:p])
        x[i, 2] <- (xb[i, 2] - betaest[, 2] %*%  x[i, 3:p])
        x[i, ] <- Vroot %*% x[i, ]
    }else if (linear == 1){
           x[i, 3:p] <- betaest %*%invbb %*% xb[i, ] +  as.matrix(D) %*% as.matrix(v^(-1/2) * (u[3:p] - m))
        x[i, 1] <- (xb[i, 1] - betaest[, 1] %*%  x[i, 3:p])
        x[i, 2] <- (xb[i, 2] - betaest[, 2] %*%  x[i, 3:p])

    }else{
        x[i, 3:p] <- (betaest %*% invbb %*%  (xb[i, ]^2 - diag(tbb))) +  as.matrix(D) %*% as.matrix(v^(-1/2) * (u[3:p] - m))
        x[i, 1] <- (xb[i, 1] - betaest[, 1] %*%  x[i, 3:p])
        x[i, 2] <- (xb[i, 2] - betaest[, 2] %*%  x[i, 3:p])
    }
}
   #     u <- runif(p, -sqrt(3), sqrt(3) )
    #    bu <- sigma %*% u
    #x[i, ] <- beta * as.numeric(invbb) * xb[i] + bu

   y <- rnorm(n, sin(2  * xb[, 1]), sqrt(d * log( 2 + xb[, 2] ^2)))
    mxy <- cbind(y, x)
    return(mxy)
}

#- (y- sin(2 * x) )^2/(2 * log(2 + x^2)) + log(1 / sqrt(2 * pi * log(2 + x^2)))
#log(e^(-(x-y)/b)/(b * (1 + e^(-(x - y)/b))^2 )) 

seff <- function(betaest, x, y, linear){
    betaest <- matrix(betaest, nrow = 3, ncol =2)
    beta <- rbind(fix, betaest)
    tbb <- t(beta) %*% beta
    tbbroot <- svd(tbb)
    tbbroot <- tbbroot$u %*% diag(tbbroot$d)^{1/2} %*% t(tbbroot$v)
    p1 <- nrow(betaest)
    p <- p1 + 2
    n <- nrow(x)
    seff <- matrix(NA, n, (p-2) * 2)
     for(i in 1: n){
        x2 = x[i, 3:p]
        xb = c(sum(x[i, ] * beta[, 1]), sum(x[i, ] * beta[, 2])) 
        invbb <- solve(t(beta) %*% beta)
        if(linear == 1){
        mx <- betaest %*% invbb %*% xb
    }else{
           mx <- betaest %*% invbb %*% (xb^2 - diag(tbb))
           mx <- c(xb[1] - betaest[, 1] %*% mx, xb[2] -  betaest[, 2] %*% mx, mx)
           mx <- (Vroot %*% mx)[3:p]
    }
        er = x2 - mx
          lf3 =c(2 * cos(2 * xb[1]) * (y[i] - sin(2 * xb[1]))/(d * log(xb[2]^2 + 2)), xb[2] * (y[i] - sin(2 * xb[1]))^2/(d * (xb[2]^2 + 2) * log(xb[2]^2 + 2)^2) - xb[2] /((xb[2]^2+ 2) * log(xb[2]^2 + 2))) #c( xb[1] * (-sin(2 * xb[1]) + y[i] - cos(2 * xb[2])) ^2 /( (xb[1] ^2 + 2) * (log(xb[1]^2 + 2))^2) + (2 * cos(2 * xb[1])) * (-sin(2 * xb[1])  + y[i] - cos(2 * xb[2]))/( log(xb[1]^2 + 2)) - xb[1] /((xb[1]^2+ 2) * log(xb[1]^2 + 2)), 2 * (cos(2 * xb[2]) - y[i] + sin(2 * xb[1])) * sin(2 * xb[2])/log(xb[1]^2 + 1)) #xb * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i]) ^2 /((xb ^2 + 2) * (log(xb^2 + 2))^2) - (-2 * cos(2 * xb) - 2 * exp(xb + 2)) * (-sin(2 * xb) - 2 * exp(xb + 2) + y[i])/(log(xb^2 + 2)) - xb /((xb^2+ 2) * log(xb^2 + 2))
        seff[i, ] <- as.vector(er %*% t(lf3))
        
    }
     mseff<- sum(apply(seff, 2, sum)^2)
    return(mseff)
}

seff3 <- function(betaest, x, y,  nonconst, t0, t1){
    
    betaest <- matrix(betaest, nrow = 3, ncol =2)
    beta <- rbind(fix, betaest)
    a11 <- t(beta[, 1]) %*% beta[, 1]
    a12 <- t(beta[, 1]) %*% beta[, 2]
    a22 <- t(beta[, 2]) %*% beta[, 2]
    b <- sqrt(3 * (a22 - a12 * a11^{-1} * a12) /pi^2) 
    k0 = a11/t0^2
    kt = k0 * t0
   
   p1 <- nrow(betaest)
    p <- p1 + 2
    n <- nrow(x)
     seff3 <- matrix(NA, n, (p-2) * 2 )
    for(i in 1: n){
        x2 = x[i, 3:p]
        xb =t(beta) %*% x[i, ] #c(sum(x[i, ] * beta[, 1]), sum(x[i, ] * beta[, 2])) 
        invbb <- solve(t(beta) %*% beta)
        mx <- betaest %*% invbb %*% xb
        er = x2 - mx
        if(nonconst == 1){
        csigm =  as.numeric(xb[1] ^2/ (t(beta[, 1]) %*% beta[, 1])) *  (diag(1, p - 2, p - 2) - betaest %*%invbb %*% t(betaest))
         sigm <-  csigm
         invcsigm = ginv(csigm)
          invsigm <-  invcsigm
        dsigma <- as.numeric(2 * xb[1] /(t(beta[, 1]) %*% beta[, 1]))* (diag(1, p - 2, p - 2) - betaest %*%invbb %*% t(betaest))
    }else{
         csigm =   (diag(1, p - 2, p - 2) - betaest %*%invbb %*% t(betaest))
         sigm <-  csigm
         invcsigm = ginv(csigm)
          invsigm <- invcsigm
    }
        dmdxb = betaest %*% invbb
       # lf1 = (xb[1] > -kt) * c(((k0-1)/(xb[1] + kt) - 1/t0) + (xb[2] - a12 * a11^{-1} * xb[1]) *  a12 * a11^{-1}/(a22 - a12^2 * a11^{-1}), - (xb[2] - a12 * a11^{-1} * xb[1]) /(a22 - a12^2 * a11^{-1}))
         lf1 = (xb[1] > -kt) * c(((k0-1)/(xb[1] + kt) - 1/t0) + (exp(xb[2]/b) - exp(a12 * a11^{-1} * xb[1]/b)) * (a12 * a11^{-1})/ (b * (exp(a12 * a11^{-1} * xb[1]/b) + exp(xb[2]/b))), -(exp(xb[2]/b) - exp(a12 * a11^{-1} * xb[1]/b))/ (b * (exp(a12 * a11^{-1} * xb[1]/b) + exp(xb[2]/b))))
     # -std ^(-2) * invbb * xb
  lf3 =c(2 * cos(2 * xb[1]) * (y[i] - sin(2 * xb[1]))/(d * log(xb[2]^2 + 2)), xb[2] * (y[i] - sin(2 * xb[1]))^2/(d * (xb[2]^2 + 2) * log(xb[2]^2 + 2)^2) - xb[2] /((xb[2]^2+ 2) * log(xb[2]^2 + 2)))
        ldmdb2 <- list()
        for(ip1 in 1:2){
        for (ip in 1:(p - 2)){
            v1 <- matrix(0, p, 2)
            v2 <- matrix(0, p-2, 2)
            v1[2 + ip, ip1] <- 1
            v2[ip, ip1] <- 1
            
            ldmdb2[[ip1 * (p - 2) + ip]] <- - t(x[i, ]) %*% (beta) %*% invbb %*% (t(v1) %*% beta + t(beta) %*% v1)%*% invbb %*% t(betaest) +  t(x[i, ]) %*% (beta) %*% invbb %*% t(v2)
        }
    }
       dmdb2 = do.call(rbind, ldmdb2)#as.numeric(xb * invbb) * diag(1, p1, p1) - 2 * as.numeric(xb * invbb^2 )* betaest %*% t(betaest)# + x2 %*% invbb %*% t(betaest) (1 - betaest^2)/(1 + betaest^2)^2 * xb
        if(nonconst == 1){
        seff3[i, ] <- as.vector(er %*% t(lf1) + mx %*% t(er) %*% invsigm %*% dmdxb + er %*% t(lf3)) + dmdb2 %*% invsigm %*% er +  as.vector(dsigma   %*% invsigm %*% er)
    }else{
         seff3[i, ] <-as.vector(er %*% t(lf1) + mx %*% t(er) %*% invsigm %*% dmdxb + er %*% t(lf3)) + dmdb2 %*% invsigm %*% er
     }
       
    }
    mseff<- sum(apply(seff3, 2, sum)^2)
    return(mseff)
    
    
}


seff42 <- function(betaest, x, y, k, t0, t1){
    
    m <-t1  * gamma(1 + 1/k)#k  * t1# -digamma(1)
    v <-t1^2* (gamma(1 + 2/k) - (gamma(1 + 1/k))^2)#k  * t1^2# pi^2 /6 
                
    betaest <- matrix(betaest, nrow = 3, ncol =2)
    beta <- rbind(fix, betaest)
    a11 <- t(beta[, 1]) %*% beta[, 1]
    a12 <- t(beta[, 1]) %*% beta[, 2]
    a22 <- t(beta[, 2]) %*% beta[, 2]
    k0 = a11/t0^2
    kt = k0 * t0
   b <- sqrt(3 * (a22 - a12 * a11^{-1} * a12) /pi^2) 
   p1 <- nrow(betaest)
    p <- p1 + 2
    n <- nrow(x)
    seff4 <- matrix(NA, n, (p-2) * 2 )
    tbb <- t(beta) %*% beta
     invbb <- solve(tbb)
     csigm =   (diag(1, p -2, p-2) - betaest %*%invbb %*% t(betaest))
    svdcsigm = svd(csigm)
    
 
    sigm <- csigm
      
    D =  svdcsigm$u %*% diag(svdcsigm$d)^(1/2) %*% t(svdcsigm$v)

    invcsigm = ginv(csigm)
    invsigm <-  invcsigm
    dmdxb = betaest %*% invbb
    invD <- ginv(D)
    # svdsigm <- svd(invsigm)
    #invD = svdsigm$u %*% diag(svdsigm$d)^(1/2) %*% t(svdsigm$v)
    # D <- solve(invD) #svdsigm$u %*% diag(svdsigm$d)^(-1/2) %*% t(svdsigm$v)
       mC <- matrix(NA,(p-2)^2, (p -2) * 2 )
         
    secm <- 1# mm3[1]#1#
    thirdm <- (gamma(1 + 3/k) * t1^3 - 3 * m * v - m^3)/v^(3/2)#mm3[2]#1.14#mm3[2]
    gamma1 <- gamma(1+ 1/k)
     gamma2 <- gamma(1+ 2/k)
     gamma3 <- gamma(1+ 3/k)
     gamma4 <- gamma(1+ 4/k)
    fourm <- (-6 * gamma1^4 + 12 * gamma1^2 * gamma2- 3 * gamma2^2 - 4 * gamma1 * gamma3 + gamma4)/(gamma2 - gamma1^2)^2 + 3#mm3[3]#12/5 + 3#mm3[3]
    iden <- diag(1, p-2, p-2)
   
    
   
    eem <- diag(secm, p-2, p-2)
    
    efv2 <- efv1 <- matrix(0, p-2, p-2)
    efv1[1, 2] <- -1
    efv1 [2, 1] <- -1
    efv2[1, 3] <- -1
    efv2[3, 1] <- -1
    
    nvI<- as.vector(diag(-1, p-2, p-2))
    efI <- matrix(0, (p - 2)^2, (p-2)^2)
    efI[, c(1, 5, 9)]<- nvI
    efe <- matrix(0, (p-2)^2, p-2)
    efe[1, 1] <- efe[5, 2] <- efe[9, 3] <- -1

    v1em <- vector('list')
    efv <- v1v1m <- vector('list')
    for(j in 1: (p -2)){
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
        v1em[[j]] <- matrix(0, (p - 2), (p - 2))
        v1em[[j]][j, j] <- thirdm
        v1v1m[[j]] <- matrix(0, (p - 2), (p - 2))
        v1v1m[[j]][j, j] <- fourm
        
        diag(v1v1m[[j]]) <- (secm) ^2
         v1v1m[[j]][j, j] <- fourm
    
        
    }
    v1em <- do.call(rbind, v1em) #(p-1)^2 by p-1
    v1v1m <- as.matrix(bdiag(v1v1m))
    Evv <- v1v1m - as.vector(iden) %*% t(as.vector(iden))-   v1em %*% t(v1em)
    
    invEvv <- ginv(Evv)
    
    efv <- do.call(rbind, efv)


for(j in 1:2){
   for(i in 1: (p - 2)){
        temp <- temp0 <- betaest
        temp[i, j] <-  betaest[i, j] + eps
        betatemp <- rbind(fix, temp)
         tbbtemp <- t(betatemp) %*% betatemp
        invbbtemp <- ginv(tbbtemp)
          csigm1 =   (diag(1, p -2, p-2) - temp %*%invbbtemp %*% t(temp))
         svdcsigm1 = svd(csigm1)
       
        D1 = std * svdcsigm1$u %*% diag(svdcsigm1$d)^(1/2) %*% t(svdcsigm1$v)
   
    
        invD1 = ginv(D1)
        DD = (invD1 - invD)/eps
 
       
        uppD<- t(D) %*% DD
        mC[, (j-1) * (p-2) + i]  <- as.vector((uppD))
        
    
    }

    
    k4 = t(mC) - t(dmdxb) %*% t(invD) %x% D
    k5 <- -k4 %*% (efv - efI - efe %*% t(v1em)) %*% invEvv
}
    for(i in 1: n){
    xb = t(beta) %*% x[i, ]
  # (p-1)^2 by (p-1)^2
   
    mx <-  betaest %*%invbb %*% xb
    er = x[i, 3:p] - mx
   # lf1 = (xb > -kt) *((k0-1)/(xb + kt) - 1/t0)
  #  lf1 = - (xb/b) - 2 * log(1 + exp(-xb/b))
   lf1 = (xb[1] > -kt) * c(((k0-1)/(xb[1] + kt) - 1/t0) + (exp(xb[2]/b) - exp(a12 * a11^{-1} * xb[1]/b)) * (a12 * a11^{-1})/ (b * (exp(a12 * a11^{-1} * xb[1]/b) + exp(xb[2]/b))), -(exp(xb[2]/b) - exp(a12 * a11^{-1} * xb[1]/b))/ (b * (exp(a12 * a11^{-1} * xb[1]/b) + exp(xb[2]/b))))
    #  lf1 = (xb[1] > -kt) * c(((k0-1)/(xb[1] + kt) - 1/t0) + (xb[2] - a12 * a11^{-1} * xb[1]) *  a12 * a11^{-1}/(a22 - a12^2 * a11^{-1}), - (xb[2] - a12 * a11^{-1} * xb[1]) /(a22 - a12^2 * a11^{-1}))
     # -std ^(-2) * invbb * xb
   lf3 =c(2 * cos(2 * xb[1]) * (y[i] - sin(2 * xb[1]))/(d * log(xb[2]^2 + 2)), xb[2] * (y[i] - sin(2 * xb[1]))^2/(d * (xb[2]^2 + 2) * log(xb[2]^2 + 2)^2) - xb[2] /((xb[2]^2+ 2) * log(xb[2]^2 + 2)))
    normer = invD %*% er
    v <- as.vector( normer %*% t(normer)  - diag(1, p - 2, p-2)) - v1em %*% normer
    
     ldmdb2 <- list()
        for(ip1 in 1:2){
        for (ip in 1:(p - 2)){
            v1 <- matrix(0, p, 2)
            v2 <- matrix(0, p-2, 2)
            v1[2 + ip, ip1] <- 1
            v2[ip, ip1] <- 1
            
            ldmdb2[[ip1 * (p - 2) + ip]] <- - t(x[i, ]) %*% (beta) %*% invbb %*% (t(v1) %*% beta + t(beta) %*% v1)%*% invbb %*% t(betaest) +  t(x[i, ]) %*% (beta) %*% invbb %*% t(v2)
        }
    }
       dmdb2 = do.call(rbind, ldmdb2)#
  
   
   k2 = - (t(dmdxb) %x%  mx + dmdb2)  %*% t(invD) #p-d by p
  #  k2 =-( mx %*% t(er)%*% invsigm %*% dmdxb + dmdb2 %*% t(invD) %*% normer)
   
  k3 <- k2 %*% t(v1em) %*% as.matrix(invEvv)
    
    seff4[i, ] <- as.vector(er %*% t(lf1) + er %*% t(lf3)) - k2   %*% normer + k3 %*% v - k5 %*% v
}
    mseff<- sum(apply(seff4, 2, sum)^2)
    return(mseff)
}





estimation <- function(beta0, betaini, n, p, itr){
    set.seed(itr + 3016)
    
    simdata <- simu2(n,  p, beta0,  k, nonconst, linear,   t0, t1)
    y <- simdata[, 1]
    x <- simdata[, 2: (p+ 1)]
    temp1 <-optim(betaini, seff, gr = NULL, x, y,linear, method = 'Nelder-Mead', control = list(abstol = 1e-6, reltol = 1e-6 ))# try(nlm( seff,betaini,  x, y))### #spg(betaini, seff, gr=NULL, method=3, lower=-Inf, upper=Inf,  project=NULL, projectArgs=NULL,  control=list(), quiet=FALSE, alertConvergence=TRUE, x, y)$par##dfsane(betaini, seff, method=2, control=list(trace = FALSE, tol= 1e-4),quiet=FALSE, alertConvergence=TRUE, x, y)$par###
    estseff <- temp1$par
   # temp1 <-try(nlm( seff,estseff,  x, y))# optim(c( estseff), seff, gr = NULL, x, y, method = 'BFGS')#
  #  estseff <- temp1$par
  #  temp1 <-optim(betaini, seff, gr = NULL, x, y, 1, method = 'BFGS')
 #   estseff1 <- temp1$par
    temp2 <-optim(c( betaini), seff3, gr = NULL, x, y,  nonconst, t0, t1,  method = 'Nelder-Mead', control = list( abstol = 1e-6, reltol = 1e-6))
#    temp2 <-optim(betaini, seff3, gr = NULL, x, y, std, std1, 0, t0, t1,   method = 'BFGS')
                                        #nlm(seff3,betaini, x, y, std, 1)# ######spg(betaini, seff3, gr=NULL, method=3, lower=-Inf, upper=Inf,  project=NULL, projectArgs=NULL,  control=list(), quiet=FALSE, alertConvergence=TRUE, x, y, std)$par#dfsane(betaini, seff3, method=2, control=list(trace = FALSE, tol = 1e-4), quiet=FALSE, alertConvergence=TRUE, x, y, std)$par#
    estseff3 <- temp2$par
    temp3 <- optim(c( betaini), seff42, gr = NULL, x, y,  k, t0, t1,  method = 'Nelder-Mead', control = list(abstol = 1e-6, reltol = 1e-6 ))
  #  temp3 <- optim(temp3$par * 1.01, seff42, gr = NULL, x, y, std, std1, k, t0, t1,   method = 'BFGS')#nlm(seff41,c(betaini), x, y, std)####optim(c(1, estseff), seff4, gr = NULL, x, y, std, method = 'L-BFGS-B')##  spg(c(1, betaini), seff4, gr=NULL, method=3, lower=-Inf, upper=Inf,  project=NULL, projectArgs=NULL,  control=list(), quiet=FALSE, alertConvergence=TRUE, x, y, std)$par#, method = 'L-BFGS-B', lower = c(1, rep(-Inf, p-1)), upper = c(1, rep(Inf, p-1)))$par
    estseff4 <- temp3$par
   # estseff4 <- estseff4[-1]/estseff4[1]
    print(c(temp1$convergence, temp2$convergence, temp3$convergence))
    print(c(round(temp1$par, 5), round(temp2$par, 5), round(temp3$par, 5)))#, temp3$convergence))
    #estseff4 <- (estseff4/estseff4[1])[-1]
    return(c(estseff, estseff3, estseff4))#,  estseff4))
}


nsim <- 1200
p <- 5
res <- matrix(NA, nsim, (p - 2) * 3 * 2)
res1 <- matrix(NA, nsim, (p-2) * 3 * 2)
res2 <- matrix(NA, nsim, (p-2) * 3 *2 )
t0 <- 0.5
t1 <- 2
#res <- rbind(res, res1)
fix <- diag(1, 2, 2)
k <- 1
simdata <- simu2(10000,  p, beta0,  k, 0, 2,   t0, t1)
x <- simdata[, 2: (p+ 1)]
normv <- cov.rob(x)$cov
Vroot <- svd(solve(normv))
Vroot <- Vroot$u %*% diag(Vroot$d)^{1/2}%*% t(Vroot$v)
beta0 <- rbind(fix, matrix(c(0.2, 0.3, 0.4, 0.3, 0.2, 0.4), nrow = 3, ncol = 2))
betaini <- as.vector(c( beta0[3:p, ] * 1.21))
    

 for(itr in 1:nsim){
    print(itr)
    std <- 1
    std1 <- sqrt(1)
    n <- 50
    
   
    temp<- try(estimation(beta0, betaini, n, p, itr))
    if(class(temp) != 'try-error')
        res[itr, ] <- temp
}
 ix <- 1000
mbeta0 <- matrix(as.vector(beta0[3:p, ]), ix, ncol = (p-2) * 2, byrow = T)
sqereff <- sqrt(apply((res[1:ix, 1:((p-2) * 2)] - mbeta0)^2, 1, mean, na.rm = T))
sqereff3 <- sqrt(apply((res[1:ix,(((p-2) * 2)+ 1):((p-2) * 2 * 2)] - mbeta0)^2, 1, mean, na.rm = T))
sqereff4 <- sqrt(apply((res[1:ix, (((p-2) * 2 * 2) + 1):((p-2) * 2 * 3)] - mbeta0)^2, 1, mean, na.rm = T))#apply((res[, 11:15] - mbeta0)^2, 1, sum)
sqereff <- sqereff[which(sqereff<= quantile(sqereff, 0.9, na.rm = T))]
sqereff3 <- sqereff3[which(sqereff3<= quantile(sqereff3, 0.9, na.rm = T))]
sqereff4 <- sqereff4[which(sqereff4<= quantile(sqereff4, 0.9, na.rm = T))]
pdf('temp.pdf')
boxplot(as.vector(sqereff), as.vector(sqereff3), sqereff4,   range = 0.8, outline = F, names = c("a", "b", 'c'), pars = (list(boxwex = 0.5)))
dev.off()

resab <- res
res3 = cbind(resab, res)
res1 <- res

abb1 = round(abs(apply(res3[1:ix, 1:3], 2, median, na.rm = T) - as.vector(beta0[3:p, ])), 4)
abb2 = round(abs(apply(res3[1:ix, 4:6], 2, median, na.rm = T) - as.vector(beta0[3:p, ])), 4)
abb3 = round(abs(apply(res3[1:ix, 7:9], 2, median, na.rm = T) - as.vector(beta0[3:p, ])), 4)

mad1 = round(apply(res3[1:ix, 1:3], 2, mad), 4)
mad2 = round(apply(res3[1:ix, 4:6], 2, mad), 4)
mad3 = round(apply(res3[1:ix, 7:9], 2, mad), 4)

print(noquote(paste(abb1, '(', mad1, ')','&',  abb2, '(', mad2, ')','&', abb3, '(', mad3, ')', '\\')))

upper.trim.mean <- function(x,trim) {
 
  mean(x[x< quantile(x, trim, na.rm = T)],  na.rm = T)
}

mvar <- function(fres){
  
    
    bias<- round(matrix(abs(apply(fres[1:1000, ], 2, median, na.rm = T) - rep(as.vector(beta0[3:p,]), 3)), nrow = (p-d) * d) , 5)
std <- round(matrix(apply(fres[1:1000, ], 2, mad, na.rm = T), nrow = (p-d) * d), 4)
    list(bias, std)
   
}
mfres <- list(fres1, fres2, fres3, fres1gamma, fres2gamma, fres3gamma)
mfres2d <- list(fres2d1, fres2d2, fres2d3)
mbiasstd <- list()
for(i in 1:3 ){
    mbiasstd[[i]] <- mvar(mfres2d[[i]])
}


print(noquote(paste(abb1, '(', mad1, ')','&',  abb2, '(', mad2, ')','&', abb3, '(', mad3, ')', '\\')))
i <- 3
print(noquote(paste(mbiasstd[[i]][[1]][, 1], ' (', mbiasstd[[i]][[2]][, 1], ')', '&',mbiasstd[[i]][[1]][, 2], ' (', mbiasstd[[i]][[2]][, 2], ')', '&', mbiasstd[[i]][[1]][, 3], ' (', mbiasstd[[i]][[2]][, 3], ')', '\\',  sep = '' )))

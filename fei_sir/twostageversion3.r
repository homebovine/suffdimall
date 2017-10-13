library(survRM2)
library(survAUC)
library(gsDesign)
library(parallel)
library(numDeriv)
library(MASS)

###use two values at me 2 * 12 and 10 *12 to get the weibull distribution that fit the data

## getweibullpara <- function(param, values, shape = NULL){
##     if(length(param) == 2){
##         shape <- param[2]
##         scale <- param[1]
## }else{
##     scale <- param[1]
## }
##     sum((1 - pweibull(c(2, 10) * 12, shape, scale ) - c(values[1],values[2]))^2)
    
## }
getweibullpara <- function(tscale, values, p, tshape){
        p1density <-  (1 - pweibull(2 * 12, scale = tscale[1], shape = tshape[1]))/(1 - pweibull(6, scale = tscale[1], shape = tshape[1]))
        p0density <- 1 - pweibull(2 * 12, scale = tscale[2], shape = tshape[2])
        d1 <- p1density * p +  p0density * (1 - p)
         p1density <- (1 - pweibull(10 * 12, scale = tscale[1], shape = tshape[1]))/(1 - pweibull(6, scale = tscale[1], shape = tshape[1]))
        p0density <- 1 - pweibull(10 * 12, scale = tscale[2], shape = tshape[2])
        d2 <- p1density * p +  p0density * (1 - p)
        sum((c(d1, d2) - values)^2)
        
    }
##get parameter for binary response
getpparam <- function(para2, para1, value, wparam){
    integrand <- function(t){
        pweibull(t - 6, para2, para1) * dweibull(t, wparam[1], wparam[2])
    }
    (integrate(integrand, 6, Inf)$value - value)^2
}


singlefulllikelihood <- function(i, data, tscale, tshape, p){
    X <- data[i, 2]
    d <- data[i, 3]
    r <- data[i, 1]
    ixcr <- data[i, 4]
 
    if(d  == 1 & ixcr == 1){
        if(r == 1){
            likelihood = dweibull(X, scale = tscale[1], shape = tshape[1])/(1 -pweibull(6, scale = tscale[1], shape = tshape[1]) ) * p
        }else{
            likelihood = dweibull(X, scale = tscale[2], shape = tshape[2]) *(1 - p)
        }
    }else if (ixcr == 1&  d  == 0){
        if(r == 1) {
            if(X <= 6){
                likelihood = 1 * p
            }else{
                likelihood =   (1 - (pweibull(X, scale = tscale[1], shape = tshape[1])  - pweibull(6, scale = tscale[1], shape = tshape[1]))/(1 -pweibull(6, scale = tscale[1], shape = tshape[1]) )) * p
            }
        }
        else{
            likelihood =  (1 -  pweibull(X, scale = tscale[2], shape = tshape[2])) * (1 - p)
        }
        
        
    }else if (ixcr == 0 & d == 1){
        likelihood =   (X >= 6) * dweibull(X, scale = tscale[1], shape = tshape[1])/(1 -pweibull(6, scale = tscale[1], shape = tshape[1]) ) * p + dweibull(X, scale = tscale[2], shape = tshape[2]) *(1 - p) 
    }else if(ixcr == 0 & d == 0){
        likelihood = min( (1 - (pweibull(X, scale = tscale[1], shape = tshape[1])  - pweibull(6, scale = tscale[1], shape = tshape[1]))/(1 -pweibull(6, scale = tscale[1], shape = tshape[1]) )), 1)   * p + (1 -  pweibull(X, scale = tscale[2], shape = tshape[2]))  * (1 - p)
        
        
    }
    if(is.nan(log(likelihood + 1e-6))){
     #  browser()
    }
    return(-log(likelihood + 1e-6))
    
    
}


## singlefulllikelihood <- function(i, data, tscale, tshape, p){
##     X <- data[i, 2]
##     d <- data[i, 3]
##     r <- data[i, 1]
##     ixcr <- data[i, 4]
 
##    if (d == 1 ){
##         likelihood = (X >= 6) *  dweibull(X, scale = tscale[1], shape = tshape[1])/(1 -pweibull(6, scale = tscale[1], shape = tshape[1]) ) * p + dweibull(X, scale = tscale[2], shape = tshape[2]) *(1 - p) 
##     }else{
##         likelihood = min( (1 - (pweibull(X, scale = tscale[1], shape = tshape[1])  - pweibull(6, scale = tscale[1], shape = tshape[1]))/(1 -pweibull(6, scale = tscale[1], shape = tshape[1]) )), 1)   * p + (1 -  pweibull(X, scale = tscale[2], shape = tshape[2]))  * (1 - p)
        
        
##     }
##     if(is.nan(log(likelihood + 1e-6))){
##      #  browser()
##     }
##     return(-log(likelihood + 1e-6))
    
    
## }


jointlikelihood <- function(tscale, data, tshape, p){
  # tscale  <- exp(tscale)
    (sapply(1:nrow(data), singlefulllikelihood, data,tscale, tshape, p)) 
}
sjointlikelihood <- function(tscale, data, tshape, p){
  # tscale <- exp(tscale)
    res <-sum(sapply(1:nrow(data), singlefulllikelihood, data, tscale, tshape, p), na.rm = T)
#print(res)
return(res)

}

wrapsjointlikelihood <- function( tpscale, data, tshape){
                                        # tscale <- exp(tscale)
    p <- tpscale[3]
    tscale <- tpscale[1:2]
   sjointlikelihood(tscale, data, tshape, p)
}

wrapjointlikelihood <- function( tpscale, data, tshape){
                                        # tscale <- exp(tscale)
    p <- tpscale[3]
    tscale <- tpscale[1:2]
   jointlikelihood(tscale, data, tshape, p)
}


pD <- function(tpscale,  X,  ind, tshape){
    tscale <-(tpscale[1:2])
    p <- tpscale[3]
  
    dervt <- (1 -  (pweibull(24, shape = tshape[1], scale = tscale[1]) - pweibull(6, shape = tshape[1], scale = tscale[1]))/(1 - pweibull(24, shape = tshape[1], scale = tscale[1]))) * p + (1 - pweibull(24, shape = tshape[2], scale = tscale[2])) * (1 - p)
    return(c(dervt))
}



#mbeta <- rep(NA, lsrv)
nsim <- 1000
simun <- 180
delta <- 0
delta1 <- 0
delta10 <- 0.05#0.05##simulation delta
delta11 <- 0.06#0.05##simulation delta
deltam1 <- c(0, 0.05, 0.05, 0.1, 0.1, 0, 0.05,   0,  0.1)
deltam2 <- c(0, 0.06, 0.12, 0.06, 0.12, 0.06, 0,   0.12, 0)
deltam <- cbind(deltam1, deltam2)
claimd0 <- 0.06# -log(1 - 0.06/0.91)##NI delta
claimd1 <- (0.05)##NI delta
msrv <- expand.grid(seq(0.82, 0.91, 0.03), seq(0.71, 0.8, 0.03))
lsrv <- nrow(deltam)
resstop <- vector("list")#matrix(NA, lsrv, 3)
resstage <- vector("list")#rep(NA, lsrv)


simufun <- function(j){
    set.seed(2020 )
    delta11 <- deltam[j, 1]
    delta10 <- deltam[j, 2]
    srvsimu <- 0.91 - delta10#the survival probability at 2 years in the simulation
    lsrvp <- 0.64- 1.5 * delta10##the survival probability at 10 years in the simulation
    srv <- 0.91 -delta##the survival probability at 2 year under the h0
    para1 <- 2##the scale parameter for the probability in the response
                                        #   Chaz0 <- (srv)n
    
    if(j  %in% c(1,2,4,  6,7, 9)){
        tshape <-c(1, 1) } else{tshape = c(1,1)}
    ## if(j == 1){
    ##     tshape = c(0.5, 0.5)
    ## }else{
    ##     tshape <- c(1, 0.9)
    ## }
                                        #c(1, 0.9)#j = 1,2, 6,7, 8, 9,  10 (1, 0.9), j = 5, 3 (0.5, 0.5)
  #  rho <- 30
  
    smalln <- 5
    mdata <- vector('list')
    nullmeanX <-0.8# #the response rate under the null
   
 
    nullmeanXsimu <-0.8 - delta11 #the response rate in the simulation
   
  
                                        #  nullmeanTsimu <- 1/nullambda
    if(j == 4|j == 5){
        itscale <- tscale <- optim(par = c(200, 10), fn = getweibullpara, gr = NULL,  values = c(srvsimu,lsrvp), p = nullmeanXsimu, tshape = tshape, method ='L-BFGS-B', lower = c(10, 10), upper = c(1000, 1000))$par #
    }else{
        itscale <- tscale <- optim(par = c(100, 10), fn = getweibullpara, gr = NULL,  values = c(srvsimu,lsrvp), p = nullmeanXsimu, tshape = tshape, method ='L-BFGS-B', lower = c(10, 10), upper = c(1000, 1000))$par
    } ##get shape and score parameter for the weibull distribution under the null (not change when simulation parameters changed)
 #   wparam <- optim(c(3), getweibullpara, gr = NULL, method ='L-BFGS-B',  lower = 0.01, upper = 10000,  values = c(srvsimu,lsrvp), pmusigma[2])$par##fix shape parameter, and obtain the scale parameter to achieve certain survival time in the simulation
 #   pparam <- optim(c(0.5), getpparam, gr = NULL, method ='L-BFGS-B', lower = 0.01, upper = 10000,  para1 = para1, value =nullmeanXsimu, c( pmusigma[2], wparam))$par ##get the scale parameter when fix the shape parameter in the binary response probability
    
  


    Kmax <- 5 * 12
    interim <- 12 * 3
    estop <- nstage <- rep(NA, nsim)
##simulation when shape parameters are all fixed
    for(itr in 1:nsim){
        data <- NULL
        for(stage in 1:interim){
            cr <- rbinom(smalln, 1, nullmeanXsimu)
            u <- runif(smalln)
            pg6 <-  pweibull(6, shape = tshape[1], scale = tscale[1])
            T <- qweibull(u * (1 - pg6) + pg6, shape = tshape[1], scale = tscale[1])
            T[cr == 0] <- rweibull(smalln, shape = tshape[2], scale = tscale[2])[cr == 0]
            
           
            respstage <- rep(stage + 6 -1, smalln)
            inistage <- rep(stage, smalln)
            C <- rexp(smalln,1/5000)
            
            data <- rbind(data, cbind(cr, T,C,  respstage ,inistage))
        }
   
        mdata[[itr]] <- data
    }

    start <-  7
    startPFS <- 3 * 12
    lowerbound <- gsDesign(Kmax - start + 1, 2,alpha = 0.1)$lower$bound
    lowerbound <- -gsDesign(Kmax - start + 1, 2,alpha = 0.05, sfupar = 1)$upper$bound
   #lowerbound <- gsDesign(Kmax - start + 1, 2,alpha = 0.1)$lower$bound##standarad
    upperbound <- gsDesign(Kmax - start + 1, 1,alpha = sqrt(0.05), sfupar = -15)$upper$bound
   
    mZ <- rep(NA, nsim)

    for(itr in 1:nsim){
       
        
        
        for(stage in start:Kmax){
            subdata  <- mdata[[itr]]
            libound <- lowerbound[stage - start + 1]
            uibound <- upperbound[stage - start + 1]
            subdata <- subdata[(stage - subdata[, 'inistage']) >= 1  , ]
            ixcr <- (stage - subdata[, 'inistage']) >= 6
          #  ixpfs <- (stage - subdata[, 'inistage']) >= 24
            ncr <- nrow(subdata)
 
            rnum <- 0
            X <- pmin(subdata[, 2], subdata[, 3], stage - subdata[, 'inistage']) + rnum
            
            ind = X == (subdata[, 2] + rnum)
           
            data <- subdata[, ]
         
            p1 <- mean(data[, 1]  & ixcr)/mean(ixcr)
           
            lower <- c(itscale[1] - 200,10, 0.4)#c(itscale[1]/30 - 150/30,10/30,0.4)
            upper <- c(1e5, 1e5,  0.99)
            optimres <- try(optim(par =c(itscale, p1), fn = wrapsjointlikelihood,  data = cbind(data[, 1] , X, ind,  ixcr), method = 'L-BFGS-B', lower = lower, upper = upper,   tshape = tshape, control = , hessian =FALSE))
            if(class(optimres) == 'try-error'){
                optimres$convergence == 1
                
            }
            
            if(class(optimres) != 'try-error'&optimres$convergence == 0 & min(abs(optimres$par -lower), abs(optimres$par -upper)) >0 ){
                
        #        hL <- try(hessian(sjointlikelihood, optimres$par, data = cbind(data[, 1] , X, ind,  ixcr), tshape = tshape, p = p1))
                hLp <- try(hessian(wrapsjointlikelihood, c(optimres$par), data = cbind(data[, 1] , X, ind,  ixcr), tshape = tshape, method.args = list(eps = 1e-4, d = 1e-6)))
                sL <- try(jacobian(wrapjointlikelihood, c(optimres$par), data = cbind(data[, 1] , X, ind,  ixcr), method = "simple",  tshape = tshape))
                sc24score<- ginv(hLp) %*% t(sL)
                XY <- sc24score %*% t(sc24score)#-sL %*% solve(hL) - ((data[, 1]  & ixcr)/mean(ixcr) - p1) %*%hLp %*% solve(hL ) 
         #   sL <- try(jacobian(jointlikelihood, optimres$par, data = cbind(data[, 1] , X, ind,  ixcr),  tshape = tshape, p = p1))
            tscale <- optimres$par
                Dss <- jacobian(pD, c(optimres$par),  X= X, ind = ind, tshape = tshape)
                Dss <- rbind(Dss, c(0, 0, 1))
                XY <- Dss %*% XY %*% t(Dss)
                est <-  min( (1 - (pweibull(2 * 12, scale = tscale[1], shape = tshape[1])  - pweibull(6, scale = tscale[1], shape = tshape[1]))/(1 -pweibull(6, scale = tscale[1], shape = tshape[1]) )), 1)   * tscale[3] + (1 -  pweibull(2 * 12, scale = tscale[2], shape = tshape[2]))  * (1 -  tscale[3]) - srv
                est <- est - XY[1, 2] * ((XY[2, 2]+  1e-16)^{-1})* (tscale[3] - p1)
                estv <- sqrt(XY[1, 1] - XY[1, 2] * ((XY[2, 2]+  1e-16)^{-1}) * XY[2, 1] )
           
            Z1 <- (tscale[3] - nullmeanX)/sqrt(XY[2, 2]) #(p1 - nullmeanX  )/ sqrt(var((data[, 1]  & ixcr)/mean(ixcr))/ncr+ 1e-6)
            Z0 <- est/estv
            if(est == 0 ){
                Z0 <- 0
            }
            if(stage <= 6){
                est <- Charz - Chaz0 
                estv <- sqrt(varYY) 
                Z <- est/estv
            }else{
                Z <- min(Z0, Z1)
            }
                #Z <- Z0 ##standard
                Z01 <- Z0 + claimd0/estv
                Z11 <- Z1 + claimd1/sqrt(XY[2, 2])#  sqrt(var((data[, 1]  & ixcr)/mean(ixcr))/ncr+ 1e-6)
                
                if(class(sc24score)!= 'try-error' & stage > 6){
                    Z11 <- min(Z01, Z11)
                }else if(stage <= 6 |is.na(optimres$par[1])){
                    Z11 <- Z01
                }else{
                    Z11 <- Z11
                }
                #Z11 <- min(Z01)#standard
                                        # + delta
        
                
            }else{
                #browser()
            #    print(c(stage, "inelse"))
                
            Z1 <- (p1 - nullmeanX  )/ sqrt(var((data[, 1]  & ixcr)/mean(ixcr))/ncr+ 1e-16)
            Z <- Z1
            Z11 <- Z1 + claimd1/  sqrt(var((data[, 1]  & ixcr)/mean(ixcr))/ncr+ 1e-16)
            Z0 <- 0
                
         
                #Z11 <- min(Z01)#standard
                                        # + delta
        
        }
            
            if((Z <= libound|stage == Kmax)& stage >= start ){
                print(c(Z0, Z1, optimres$par))
                nstage[itr] <- stage
                mZ[itr] <- Z
                if(Z <= libound){
                    estop[itr] <- 1
                }else if(Z11 > qnorm(1 - sqrt(0.15))){#qnorm(1 - (0.05))#standard
                    estop[itr] <- 2
                }else{
                    estop[itr] <- 0
                }
                break
            }
           # print(stage)
        }
        print(c(itr, stage))
    }
   return(list(estop, nstage))
    ## resstop[j, 1] <- mean(estop ==1)
    ## resstop[j, 2] <- mean(estop == 2)
    ## resstop[j, 3] <- mean(nstage < 60)
    ## resstage[j] <- mean(nstage)
}

res <- mclapply(1:9, simufun, mc.cores =detectCores(all.tests = FALSE, logical = TRUE))
save(res, deltam, file = "feidesign8-3.Rdata")
## resstop <-  matrix(NA, 9, 3)
## resstage <- ressample <- matrix(NA, 9, 5)
## for(j in 1:9){
##     estop <- res[[j]][[1]]
##     nstage <- res[[j]][[2]]
##     ix <- nstage > 36
##     nsample <- nstage * 5
##     nsample[ix] <- 36 * 5
##     resstop[j, 1] <- mean(estop ==1)
##     resstop[j, 2] <- mean(estop == 2)
##     resstop[j, 3] <- mean(nstage < 60)
##     resstage[j, ] <-c(mean(nstage), sd(nstage), quantile(nstage, c(0.25, 0.5, 0.75)))
##     ressample[j, ] <-c(mean(nsample), sd(nsample), quantile(nsample, c(0.25, 0.5, 0.75)))
## }

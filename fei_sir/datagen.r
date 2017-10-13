library("MASS")
#library("expm")
p <- 50
n <- 300
d <- 2
k <- 6
nonnorm = 0
r = 0.5^abs((1:p) %*% t(rep(1, p)) - rep(1, p)%*%  t(1:p))
r = r * p
r1 = 0.5^abs((1:k) %*% t(rep(1, k)) - rep(1, k)%*%  t(1:k))
#diag(r) = p
nsim <- 1200
resmat <- matrix(NA, nsim * n, (k))
resmat3 <- resmat4 <- resmat6 <- resmat7 <- rep(NA, nsim * n)
beta1 <- diag(c(1/sqrt(6), -1/sqrt(6)), d, d)#rep(1, k)/sqrt(k)#
#beta2 <- beta1
beta2 <- matrix(c(1/sqrt(6), 1/sqrt(6), 1/sqrt(6), 1/sqrt(6), 1/sqrt(6), -1/sqrt(6),  1/sqrt(6), -1/sqrt(6)), 4, 2)
#beta2[seq(2, k, 2)] <- -beta1[seq(2, k, 2)]
beta <- rbind(beta1, beta2)
x = mvrnorm(n, rep(0, p), r);
## x <- (x - matrix(apply(x, 2, mean), n, p, byrow = T))%*% solve(cov(x))
eigenx <- eigen(x%*%t(x))
B <-  n^(-1)*t(x) %*% (eigenx$vector[, 1:k] * sqrt(n))
set.seed(2008)
for(itr in 1: nsim){
## o <- order(svdx$d, decreasing = TRUE)
F <-  mvrnorm(n, rep(0, k), r1)#matrix(rnorm(n * k), n, k)

if(nonnorm == 1){
F[, 3] = abs(F[, 1]) + abs(F[, 2]) + abs(F[, 1])  * rnorm(n)
F[, 4] = (F[, 1] + F[, 2])^2 + abs(F[, 2]) * rnorm(n)
F[, 5] = rbinom(n, 1, exp(F[, 2])/(1 +exp(F[, 2]) ))
F[, 6] = rbinom(n, 1, pnorm(F[, 2]))
}
#B <- matrix(rnorm(p * k), p, k)
sB <- sign(B[1, ])
B<- t(apply(B, 1, '*', sB ))
cF <- cov(F )
svdcF <- svd(cF)
sovleF <- svdcF$u %*% diag((svdcF$d)^(-1/2))%*% t(svdcF$v)
F <- (F- matrix(apply(F, 2, mean), n, k, byrow = T)) %*% sovleF
#sF <- sign(F[1, ])
#F<- t(apply(F, 1, '*', sF ))
x<- t(B%*%t(F)) + matrix(rnorm(n * p, 0, sqrt(1/(n * 2))), n, p)# mvrnorm(n, rep(0, p), r/p)
#x <- (x- matrix(apply(x, 2, mean), n, p, byrow = T))

svdx <- eigen(x%*%t(x))

hatF <- svdx$vector[, 1:k] * sqrt(n)
hatB  <- n^(-1)*t(x) %*% hatF
sB <- sign(hatB[1, ])
hatF<- t(apply(hatF, 1, '*', sB ))
cF <- cov(hatF )
svdcF <- svd(cF)
sovleF <- svdcF$u %*% diag((svdcF$d)^(-1/2))%*% t(svdcF$v)
hatF <- (hatF- matrix(apply(hatF, 2, mean), n, k, byrow = T)) %*% sovleF
#sF <- sign(hatF[1, ])
#hatF<- t(apply(hatF, 1, '*', sF ))

 

y3 <- (((F %*% beta[, 1])^2 + (F %*% beta[, 2])^2)) + 0.5 * rnorm(n)
y4 <- (F %*% beta[, 1])/(0.5 + (F%*% beta[, 2] + 1.5) ^2) + 0.5 * rnorm(n)
y6 <- (((F %*% beta[, 1])^2 + 2 * abs(F %*% beta[, 2]))) + 0.1 * abs(F %*% beta[, 2]) * rnorm(n)
y7 <- (F %*% beta[, 1])^2+  abs((F %*% beta[, 2])+ 1) + 0.1* (F %*% beta[, 2])^2 *   rnorm(n)
resmat[(((itr-1)* n + 1 ):((itr)* n)), ] <- hatF
resmat3[(((itr-1)* n + 1 ):((itr)* n))] <- y3
resmat4[(((itr-1)* n + 1 ):((itr)* n))] <- y4
resmat6[(((itr-1)* n + 1 ):((itr)* n))] <- y6
resmat7[(((itr-1)* n + 1 ):((itr)* n))] <- y7
}
write.table(resmat, file = "datasrd.tab", row.names = FALSE, col.names = FALSE)
write.table(resmat3, file = "datasrd3.tab", row.names = FALSE, col.names = FALSE)
write.table(resmat4, file = "datasrd4.tab", row.names = FALSE,
col.names = FALSE)
write.table(resmat6, file = "datasrd6.tab", row.names = FALSE,
col.names = FALSE)
write.table(resmat7, file = "datasrd7.tab", row.names = FALSE,
col.names = FALSE)
write.table(beta, file = "beta.tab", row.names = FALSE, col.names = FALSE)



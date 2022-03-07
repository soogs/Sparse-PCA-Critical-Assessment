# Script to produce the results for the toy example illustration #

library(RegularizedSCA)

source("../Functions/spca_adj_lasso.R")
source("../Functions/paper1_uslpca.R")
source("../Functions/paper1_spca_adj.R")


# bmatrix function that converts an R matrix into a Latex matrix format
bmatrix <- function(x, digits=NULL, ...) {
  library(xtable)
  default_args = list(include.colnames=FALSE, only.contents=TRUE,
                      include.rownames=FALSE, hline.after=NULL, comment=FALSE,
                      print.results=FALSE)
  passed_args = list(...)
  calling_args = c(list(x=xtable(x, digits=digits)),
                   c(passed_args,
                     default_args[setdiff(names(default_args), names(passed_args))]))
  cat("\\begin{bmatrix*}[r]\n",
      do.call(print.xtable, calling_args),
      "\\end{bmatrix*}\n")
}


# the rotated bilinear-2 loadings i want in the end:
p_final <- matrix(c(1,1,0,0,1,1), ncol = 2)

# i can do svd on p
svdp <- svd(p_final)

# P = Q S B N
q <- svdp$u
# s <- diag(svdp$d) * sqrt(5)
s <- diag(svdp$d) * sqrt(4)
b <- t(svdp$v)

# round(q %*% s %*% b * 1/(sqrt(5)),5) # equal to p_final

round(q %*% s %*% b * 1/(sqrt(4)),5) # equal to p_final

# i need a bilinear-1 situation that can also be expressed with rotation

set.seed(033)
tmat <- MASS::mvrnorm(n = 5, mu = rep(0,2), Sigma = diag(2))
tmat <- scale(tmat, T, F)
tmat <- qr.Q(qr(tmat))

tmat %*% b

x <- tmat %*% s %*% t(q)

# bilinear_t <- tmat * sqrt(5)
bilinear_t <- tmat * sqrt(4)

# bilinear_p <- q %*% s / sqrt(5)
bilinear_p <- q %*% s / sqrt(4)

sum((x - bilinear_t %*% t(bilinear_p))^2)

# standardizing x
# x_std <- scaleData(x)
x_std <- scale(x, T, T)

# extracting the tmat, with bilinear-2 formulation:
svd_bilinear2 <- svd(x_std)
# bilinear_t_obs <- svd_bilinear2$u[,1:2] * sqrt(5)
bilinear_t_obs <-  svd_bilinear2$u[,1:2] * sqrt(4) 

cor(x_std, bilinear_t_obs)
# t(x_std) %*% bilinear_t_obs / 5
t(x_std) %*% bilinear_t_obs / 4
# bilinear_p_obs <- svd_bilinear2$v[,1:2] %*% diag(svd_bilinear2$d[1:2]) / sqrt(5)
# all good. 

bilinear_p_obs <- svd_bilinear2$v[,1:2] %*% diag(svd_bilinear2$d[1:2]) / sqrt(4)

# now rotating the solutions of bilinear-2
rot_t_obs <- bilinear_t_obs %*% b
# t(x_std) %*% rot_t_obs / 5
t(x_std) %*% rot_t_obs / 4
rot_p_obs <- bilinear_p_obs %*% b

# now we can work back the example, so that 
# we can already get a standardized X, from our original x = udv
u_final <- svd_bilinear2$u[,1:2]
s_final <- diag(svd_bilinear2$d[1:2])
v_final <- svd_bilinear2$v[,1:2]

# then the simple structure P that we start with would be:
# P = V S B N

# p_rot <- v_final %*% s_final %*% b / sqrt(5)
p_rot <- v_final %*% s_final %*% b

p_rot <- p_rot[,c(2,1)]
bmatrix(p_rot)

# t_rot <- u_final %*% b * sqrt(5)
t_rot <- u_final %*% b

t_rot <- t_rot[,c(2,1)]

bmatrix(t_rot)

x_final <- t_rot %*% t(p_rot)

sum((x_final - t_rot %*% t(p_rot))^2)
bmatrix(x_final)

# eigenvalue decomposition
bmatrix(t(x_final) %*% x_final)
eigs <- eigen(t(x_final) %*% x_final)
eigs$vectors <- eigs$vectors  * -1

sum((eigs$vectors[,1:2] %*% diag(eigs$values[1:2]) %*% t(eigs$vectors[,1:2]) - t(x_final) %*% x_final)^2)

bmatrix(diag(eigs$values[1:2]))
bmatrix(eigs$vectors[,1:2])

# bilinear-1 pca through svd
svd1 <- svd(x_final)
svd1$u <- svd1$u * -1
svd1$v <- svd1$v * -1

sum((svd1$u %*% diag(svd1$d) %*% t(svd1$v) - x_final)^2)

bmatrix(x_final)
bmatrix(svd1$u[,1:2])
bmatrix(diag(svd1$d[1:2]))
bmatrix(svd1$v[,1:2])

# bilinear-1 T matrix
x_final %*% svd1$v[,1:2]
svd1$u[,1:2] %*% diag(svd1$d[1:2])

sum((x_final %*% svd1$v[,1:2] %*% t(svd1$v[,1:2]) - x_final)^2)

# bilinear-2 pca through scaling the columns of bilinear-1
bmatrix(x_final)
# bmatrix(svd1$u[,1:2] * sqrt(5))
# bmatrix(svd1$v[,1:2] %*% diag(svd1$d[1:2]) * 1/sqrt(5))

bmatrix(svd1$u[,1:2] * sqrt(4))
bmatrix(svd1$v[,1:2] %*% diag(svd1$d[1:2]) * 1/sqrt(4))

# (t(x_final) %*% svd1$u[,1:2] * sqrt(5)) / 5
# cor(x_final, svd1$u[,1:2] * sqrt(5))

(t(x_final) %*% svd1$u[,1:2] * sqrt(4)) / 4
cor(x_final, svd1$u[,1:2] * sqrt(4))

# expressing bilinear-2 by x = xwp
# w_bilinear2 <- MASS::ginv(x_final) %*% (svd1$u[,1:2] * sqrt(5))
w_bilinear2 <- MASS::ginv(x_final) %*% (svd1$u[,1:2] * sqrt(4))

bmatrix(w_bilinear2)

sum((x_final %*% w_bilinear2 %*% t(svd1$v[,1:2] %*% diag(svd1$d[1:2]) * 1/sqrt(5)) - x_final)^2)

sum((x_final %*% w_bilinear2 %*% t(svd1$v[,1:2] %*% diag(svd1$d[1:2]) * 1/sqrt(4)) - x_final)^2)

# expressign billinear-2 rotation by x = xwp
w_rot <- MASS::ginv(x_final) %*% t_rot

sum((x_final %*% w_rot %*% t(p_rot) - x_final)^2)

bmatrix(w_rot)

new_rotation <- b * -1

# p_rot %*% t(new_rotation) %*% solve(diag(svd1$d[1:2])) * sqrt(5)

p_rot %*% t(new_rotation) %*% solve(diag(svd1$d[1:2])) * sqrt(4)
# w_rot %*% t(new_rotation) %*% diag(svd1$d[1:2]) * 1/sqrt(5)
w_rot %*% t(new_rotation) %*% diag(svd1$d[1:2]) * 1/sqrt(4)

####

set.seed(113)
psparse <- spca_P_cardinality(X = x_final, R = 2, P = NULL, n_zeros = 1, MAXITER = 10000, inits = "multistart", nrstart = 100, stop_value = 1e-100, seed = 113)

min(psparse$LOSS)

psparse$bunch[[which.min(psparse$LOSS)]]

psparse$bunch[[100]]

psparse_svd_w <- MASS::ginv(x_final) %*% (psparse$bunch[[100]]$Tmat * -1)

bmatrix(psparse$bunch[[100]]$Tmat * -1)
bmatrix(psparse$bunch[[100]]$P * -1)
bmatrix(psparse_svd_w)

# vaf calculation #
12 - sum((x_final - psparse$bunch[[100]]$Tmat[,1] %*% t(psparse$bunch[[100]]$P[,1]))^2)
12 - sum((x_final - psparse$bunch[[100]]$Tmat[,2] %*% t(psparse$bunch[[100]]$P[,2]))^2)

12 - sum((x_final - psparse$bunch[[100]]$Tmat %*% t(psparse$bunch[[100]]$P))^2)

sum((x_final - psparse$bunch[[100]]$Tmat %*% t(psparse$bunch[[100]]$P))^2)

# psparse multi #
bmatrix(psparse$bunch[[which.min(psparse$LOSS)]]$Tmat[,c(2,1)])
bmatrix(psparse$bunch[[which.min(psparse$LOSS)]]$P[,c(2,1)])

psparse_multi_w <- MASS::ginv(x_final) %*% psparse$bunch[[which.min(psparse$LOSS)]]$Tmat[,c(2,1)]

bmatrix(psparse_multi_w)

# vaf calculation #

t(psparse$bunch[[which.min(psparse$LOSS)]]$Tmat) %*% psparse$bunch[[which.min(psparse$LOSS)]]$Tmat
# the components are orthogonal (uncorrelated)

sum((x_final - psparse$bunch[[which.min(psparse$LOSS)]]$Tmat %*% t(psparse$bunch[[which.min(psparse$LOSS)]]$P))^2)

sum((x_final - psparse$bunch[[which.min(psparse$LOSS)]]$Tmat[,1] %*% t(psparse$bunch[[which.min(psparse$LOSS)]]$P[,1]))^2)


sum((x_final - psparse$bunch[[which.min(psparse$LOSS)]]$Tmat[,2] %*% t(psparse$bunch[[which.min(psparse$LOSS)]]$P[,2]))^2)


psparse$bunch[[which.min(psparse$LOSS)]]$Tmat %*% t(psparse$bunch[[which.min(psparse$LOSS)]]$P)

TuckerCoef(psparse$bunch[[which.min(psparse$LOSS)]]$P, psparse$bunch[[100]]$P)


bmatrix(psparse$bunch[[100]]$Tmat %*% t(psparse$bunch[[100]]$P))

set.seed(113)

wsparse <- spca_adj(x = x_final, K = 2, para = c(2,2), type = "predictor", sparse = "varnum", inits = "multistart", nrstart = 100, eps.conv = 1e-30)

wsparse$bunch[[which.min(wsparse$LOSS)]]$Wraw*-1
bmatrix(wsparse$bunch[[which.min(wsparse$LOSS)]]$Wraw * -1)

bmatrix(wsparse$bunch[[which.min(wsparse$LOSS)]]$Pmat * -1)

wsparse_multi_t <- x_final %*% (wsparse$bunch[[which.min(wsparse$LOSS)]]$Wraw * -1)

bmatrix(wsparse_multi_t)

# for w-sparse multi, the components are correlated. So it is not really possible to calculat the vaf per component..
t(wsparse_multi_t) %*% wsparse_multi_t 

sum((x_final - x_final %*% wsparse$bunch[[which.min(wsparse$LOSS)]]$Wraw %*%  t(wsparse$bunch[[which.min(wsparse$LOSS)]]$Pmat))^2)

sum((x_final - x_final %*% wsparse$bunch[[which.min(wsparse$LOSS)]]$Wraw[,1] %*%  t(wsparse$bunch[[which.min(wsparse$LOSS)]]$Pmat[,1]))^2)

sum((x_final - x_final %*% wsparse$bunch[[which.min(wsparse$LOSS)]]$Wraw[,2] %*%  t(wsparse$bunch[[which.min(wsparse$LOSS)]]$Pmat[,2]))^2)

TuckerCoef(wsparse$bunch[[which.min(wsparse$LOSS)]]$Wraw, wsparse$bunch[[100]]$Wraw)
TuckerCoef(wsparse$bunch[[which.min(wsparse$LOSS)]]$Pmat, wsparse$bunch[[100]]$Pmat)

wsparse$bunch[[100]]$Wraw

wsparse$bunch[[100]]$Wraw <- wsparse$bunch[[100]]$Wraw * -1

wsparse$bunch[[100]]$Pmat <- wsparse$bunch[[100]]$Pmat * -1

bmatrix(wsparse$bunch[[100]]$Wraw)
bmatrix(wsparse$bunch[[100]]$Pmat)

wsparse_svd_t <- x_final %*% wsparse$bunch[[100]]$Wraw

bmatrix(wsparse_svd_t)

# with w-sparse SVD, the components are uncorrelated.
t(wsparse_svd_t) %*% wsparse_svd_t 

sum((x_final - x_final %*% wsparse$bunch[[100]]$Wraw %*%  t(wsparse$bunch[[100]]$Pmat))^2)
sum((x_final - x_final %*% wsparse$bunch[[100]]$Wraw[,1] %*%  t(wsparse$bunch[[100]]$Pmat[,1]))^2)
sum((x_final - x_final %*% wsparse$bunch[[100]]$Wraw[,2] %*%  t(wsparse$bunch[[100]]$Pmat[,2]))^2)

bmatrix(x_final %*% wsparse$bunch[[100]]$Wraw %*%  t(wsparse$bunch[[100]]$Pmat))

sum((x_final - x_final %*% wsparse$bunch[[100]]$Wraw %*%  t(wsparse$bunch[[100]]$Pmat))^2)

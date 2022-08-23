
#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestion Moderna de portafolio
# Capitulo 14: Paridad de Riesgo
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------
#
# Notas:
# 1. Requiere cargar primero la funcion para importar precios  
#    de Yahoo Finance: "func_precios"
#-----------------------------------------------------------
#-----------------------------------------------------------


rm(list=ls())
library(zoo)
library(xts)
library(quantmod)
library(quadprog)


## Acciones seleccionadas
activos <- c("ADBE","MCD","MSCI","MSFT","NEE","PG","RSG","WMT")
fechai <- '2015-12-01'
fechaf <- '2020-12-31'
periodicidad <- "monthly"   
precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]

(mu <- colMeans(retornos)*12)
(cov <- cov(retornos)*12)
var <- diag(cov)
(sigma <- sqrt(var))
n <- length(mu)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

# Ejemplo 14.2
# Paridad de Riesgo Naive

unos <- rep(1, n)
wrpnaive <- 1/sigma/(sum(1/sigma))

# Contribución al riesgo (CR):
(crnaive <- round(t(wrpnaive)*sigma,4))

# BarPlot

windows()
barplot(t(wrpnaive),beside=TRUE,col="black",ylim=c(0,0.16),ylab="Part. (%)") #

windows()
barplot(crnaive,beside=TRUE,ylim=c(0,0.025),col="black",ylab="Contribución (%)") 


# --------------------------------------------------

## Portafolio Optimo de Markowitz: PMVG

cov_inv_1 <- solve(cov, unos) 
(wpmvg <- (1/as.numeric(unos %*% cov_inv_1)) * cov_inv_1)
(crmv <- round(t(wpmvg)*sigma,4))

# Contribución relativa al riesgo (CRR):

(crrnaive <- crnaive/sum(crnaive))
(crrmv <- crmv/sum(crmv))

crcomp <- rbind(crrnaive,crrmv)

# BarPLots
windows()
barplot(crcomp,beside=TRUE,ylab="Contribución (%)",ylim=c(0,0.25))
legend("topright",legend = c("NRP","PMVG"), fill=c("black","gray"))

wcomp <- rbind(wrpnaive,wpmvg)

windows()
barplot(wcomp,beside=TRUE,ylab="Part (%)",ylim=c(0,0.25))
legend("topright",legend = c("NRP","PMVG"), fill=c("black","gray"))


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

# Ejemplo 14.3
# Risk Parity Vanilla

# Forma 1: solucion del polinomio

library(rootSolve)

b <- rep(1/n,n)
f_root <- function(x,Sigma){
    cov <- Sigma
    n <- ncol(cov)
    return(cov%*%x - b/x)
}

x_root <- multiroot(f_root, start=b, parms=cov)$root
wRPv <- x_root/sum(x_root)
names(wRPv) <- activos

# Contribución al riesgo (CR vanilla):
(crvanilla <- as.vector(wRPv * (cov%*%wRPv)))
names(crvanilla) <- activos
(crrvanilla <- crvanilla/sum(crvanilla))

# BarPlot

windows()
barplot(t(wRPv),beside=TRUE,col="black",ylim=c(0,0.16),ylab="Part. (%)") #

windows()
barplot(crrvanilla,beside=TRUE,col="black",ylim=c(0,0.14),ylab="Contribución (%)") 

# Comparación de resultados

crcomparac <- rbind(crvanilla,crnaive)

windows()
barplot(crcomparac,beside=TRUE,ylab="Contribución (%)",ylim=c(0,0.16))
legend("topright",legend = c("VRP","NRP"), fill=c("black","gray"))

wRPcomp <- rbind(wRPv,wrpnaive)

windows()
barplot(wRPcomp,beside=TRUE,ylab="Part (%)",ylim=c(0,0.2))
legend("topright",legend = c("VRP","NRP"), fill=c("black","gray"))


# Ejemplo 14.4
# Forma 2: Metodo de Newton

x0 <- rep(1/n,n)
Sigma <- cov
fn_convex <- function(x, Sigma){
    N <- nrow(Sigma)
    return(0.5 * t(x) %*% Sigma %*% x - (1/N)*sum(log(x)))
}

result <- optim(par = x0, fn = fn_convex, Sigma = Sigma, method = "BFGS")
x_cvx <- result$par
wRPv2 <- x_cvx/sum(x_cvx)
names(wRPv2) <- activos

# Contribución al riesgo (CR vanilla):
(crvanilla2 <- as.vector(wRPv2 * (cov%*%wRPv2)))
names(crvanilla2) <- activos
(crrvanilla <- crvanilla2/sum(crvanilla2))

# Comparación
wRPcomp2 <- rbind(wRPv2,wRPv)

windows()
barplot(wRPcomp2,beside=TRUE,ylab="Part (%)",ylim=c(0,0.2))
legend("topright",legend = c("VRP2","VRP1"), fill=c("black","gray"))


# --------------------------------------------------------------------
# --------------------------------------------------------------------

# Evaluacion de desempeno 

# Indice
indice <- c("^GSPC")
p.indice <- precios(indice,fechai,fechaf,periodicidad)
r.indice <- diff(log(p.indice))[-1,]
(mu.indice <- mean(r.indice))
(sigma.indice <- sd(r.indice))
(sharpe.ind <- mu.indice/sigma.indice)

# In-sample
valor <- 100 
t <- nrow(retornos)
n <- ncol(retornos)

# Retornos historicos
# PMVG
rpmv <- retornos%*%wpmvg
(mu.rpmvg <- mean(rpmv))
(sigma.pmvg <- sd(rpmv))
(sharpe.pmvg <- mu.rpmvg/sigma.pmvg)

# Paridad de Riesgo
rppr <- retornos%*%wRPv
(mu.rppr <- mean(rppr))
(sigma.pr <- sd(rppr))
(sharpe.pr <- mu.rppr/sigma.pr)

# Valor del portafolio
# PMVG
port.mv <- matrix(0, nrow=t)
port.mv[1] <- valor
for(i in 2:t){
    port.mv[i] <- port.mv[i-1]*exp(rpmv[i-1])
}

# Paridad de Riesgo
port.pr <- matrix(0, nrow=t)
port.pr[1] <- valor
for(i in 2:t){
    port.pr[i] <- port.pr[i-1]*exp(rppr[i-1])
}

# Benchmark
v.benchmark <- matrix(0, nrow=t)
v.benchmark[1] <- valor
for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.indice[i-1])
}

Performance <- cbind(port.pr,port.mv,v.benchmark)
Performance <- ts(Performance,start=2016, frequency=12)
colnames(Performance) <- c("PR","PMVG","Benchmark")

windows()
plot(Performance[,1], type='l',main="Evaluacion de desempeño", 
     ylab="Valor del portafolio", xlab="Tiempo (años)",bty="L")
lines(Performance[,2], lty=2)
lines(Performance[,3], col='darkgray')
legend("topleft",c("PR","PMVG","S&P 500"),
       lty =c(1,2,1), col=c("black","black","darkgray"))

#------------------------------------------------------------------
# Desempeno Out-sample

fechaifm <- '2020-12-01'
fechaffm <- '2021-12-31'

preciosfm <- precios(activos,fechaifm,fechaffm,periodicidad)
retornosfm <- diff(log(preciosfm))[-1,]
indicefm <- precios(indice,fechaifm,fechaffm,periodicidad)
r.indicefm <- diff(log(indicefm))[-1,]
(mu.indicefm <- mean(r.indicefm))
(sigma.indicefm <- sd(r.indicefm))
(sharpe.indicefm <- mu.indicefm/sigma.indicefm)

valor <- 100 
t <- nrow(retornosfm)+1
n <- ncol(retornosfm)

# Retornos 
# PMVG
rpmvfm <- retornosfm%*%wpmvg
(mu.rpmvgfm <- mean(rpmvfm))
(sigma.pmvgfm <- sd(rpmvfm))
(sharpe.pmvgfm <- mu.rpmvgfm/sigma.pmvgfm)

# PR
rpprfm <- retornosfm%*%wRPv
(mu.rpprfm <- mean(rpprfm))
(sigma.prfm <- sd(rpprfm))
(sharpe.prfm <- mu.rpprfm/sigma.prfm)

# Valor del portafolio
# PMVG
port.mvfm <- matrix(0, nrow=t)
port.mvfm[1] <- valor
for(i in 2:t){
    port.mvfm[i] <- port.mvfm[i-1]*exp(rpmvfm[i-1])
}

# Paridad de Riesgo
port.prfm <- matrix(0, nrow=t)
port.prfm[1] <- valor
for(i in 2:t){
    port.prfm[i] <- port.prfm[i-1]*exp(rpprfm[i-1])
}

# Benchmark
v.benchmarkfm <- matrix(0, nrow=t)
v.benchmarkfm[1] <- valor
for(i in 2:t){
    v.benchmarkfm[i] <- v.benchmarkfm[i-1]*exp(r.indicefm[i-1])
}

Performancefm <- cbind(port.prfm,port.mvfm,v.benchmarkfm)
Performancefm <- ts(Performancefm,start=c(2021,1), frequency=12)
colnames(Performancefm) <- c("PR","PMVG","Benchmark")

windows()
plot(Performancefm[,1], type='l',main="Evaluacion de desempeno", bty="L",
     ylab="Valor del portafolio", xlab="Tiempo (meses)", ylim=c(min(Performancefm), max(Performancefm)))
lines(Performancefm[,2], lty=2)
lines(Performancefm[,3], col='darkgray')
legend("topleft",c("PMVG","S&P 500"),
       lty =c(1,2,1), col=c("black","black","darkgray"))


#------------------------------------------------------------------
#------------------------------------------------------------------

## Formulación general: metodo de Feng y Palomar (2015):

Sigma <- cov
N <- n
compute_gA <- function(w, Sigma) {
    N <- length(w)
    g <- rep(NA, N^2)
    A <- matrix(NA, N^2, N)
    for (i in 1:N) {
        Mi <- matrix(0, N, N)
        Mi[i, ] <- Sigma[i, ]
        for (j in 1:N) {
            Mj <- matrix(0, N, N)
            Mj[j, ] <- Sigma[j, ]
            #g[i + (j-1)*N]   <- t(w) %*% (Mi - Mj) %*% w
            g[i + (j-1)*N]   <- w[i]*(Sigma[i, ] %*% w) - w[j]*(Sigma[j, ] %*% w)  # faster
            A[i + (j-1)*N, ] <- (Mi + t(Mi) - Mj - t(Mj)) %*% w
        }
    }
    # # g can be computed much more efficiently with this code:
    # wSw <- w * (Sigma %*% w)
    # g <- rep(wSw, times = N) - rep(wSw, each = N)
    return(list(g = g, A = A))
}

#  loop of the SCA algorithm:

library(quadprog)  # install.packages("quadprog")

# parameters
max_iter <- 40
tau <- 1e-6
zeta <- 0.1
gamma <- 0.99
# initial point
obj_value <- NULL
w_SCA <- rep(1/N, N)  # initial point

for (k in 1:max_iter) {
    # compute parameters for QP
    gA <- compute_gA(w_SCA, Sigma)
    g <- gA$g
    A <- gA$A
    Q <- 2 * t(A) %*% A + tau*diag(N)  # faster code is: crossprod(A) = t(A) %*% A
    q <- 2 * t(A) %*% g - Q %*% w_SCA
    obj_value <- c(obj_value, sum(g^2))
    
    # # solve the inner QP with CVXR
    # w_ <- Variable(N)
    # prob <- Problem(Minimize(0.5*quad_form(w_, Q) + t(q) %*% w_),
    #                 constraints = list(sum(w_) == 1))
    # result <- solve(prob)
    # w_ <- as.vector(result$getValue(w_))
    
    # solve the problem with solve.QP() which is much faster than CVXR)
    w_ <- solve.QP(Q, -q, matrix(1, N, 1), 1, meq = 1)$solution
    
    # next w
    gamma <- gamma*(1 - zeta*gamma)
    w_SCA_prev <- w_SCA
    w_SCA <- w_SCA + gamma*(w_ - w_SCA)
    
    # stopping criterion
    if (max(abs(w_SCA - w_SCA_prev)) <= 1e-4*max(abs(w_SCA_prev)))
        break
    # if (k>1 && abs(obj_value[k] - obj_value[k-1]) <= 1e-2*obj_value[k-1])
    #   break
}
cat("Number of iterations:", k)

names(w_SCA) <- activos
windows()
barplot(t(w_SCA),xlab ="stocks", ylab="Part (%)", beside = TRUE)

# Comparación de resultados
wRPcomp3 <- rbind(w_SCA,wRPv2)

windows()
barplot(wRPcomp3,beside=TRUE,ylab="Part (%)",ylim=c(0,0.2))
legend("topright",legend = c("SCA","VRP"), fill=c("black","gray"))

# ---------------------------------------------------------------------




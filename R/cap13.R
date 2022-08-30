#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestión Moderna de portafolio
# Capitulo 13: Optimización Robusta
# Copyright 2022
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

mu <- colMeans(retornos)
cov <- cov(retornos)
var <- diag(cov)
sigma <- sqrt(var)
corr <- cor(retornos)
rf <- 0 
n <- length(mu)

# Indice
indice <- c("^GSPC")
p.indice <- precios(indice,fechai,fechaf,periodicidad)
r.indice <- diff(log(p.indice))[-1,]
mu.indice <- mean(r.indice)
sigma.indice <- sd(r.indice)

##-----------------------------------------------------------
##-----------------------------------------------------------

# Ejemplo 13.1
# Solución modelo MV
# Formas convexas para PMVG y PT
library(matlab)
library(MASS)
library(CVXR)

portfolioMarkowitz <- function(mu, cov){
    w <- Variable(n)
    prob <- Problem(Minimize( quad_form(w, cov) ),
                    constraints = list(w >= 0,  ones(1,n)%*%w==1 )) 
    result <- solve(prob)
    return(as.vector(result$getValue(w)))
}
w_Markowitz <- round(portfolioMarkowitz(mu, cov),4)
names(w_Markowitz) <- colnames(cov)
w_Markowitz <- cbind(w_Markowitz)

# MV - Sharpe
portfolioMaxSharpeRatio <- function(mu, cov) {
    w <- Variable(n)
    prob <- Problem(Minimize(quad_form(w, cov)),
                    constraints = list(w >= 0, t(mu) %*% w == 1))
    result <- solve(prob)
    w <- as.vector(result$getValue(w)/sum(result$getValue(w)))
    names(w) <- colnames(cov)
    return(w)
}
w_MaxSR <- round(portfolioMaxSharpeRatio(mu, cov),4)
names(w_MaxSR) <- colnames(cov)
w_MaxSR  <- cbind(w_MaxSR)

# Resultados 
windows()
barplot(t(w_Markowitz), main="PMVG",beside = TRUE,ylab="Part (%)",
        ylim=c(0,0.3))
windows()
barplot(t(w_MaxSR), main="PMVG",beside = TRUE,ylab="Part (%)",
        ylim=c(0,0.4))

# Pesos optimos
wMV <- cbind(w_Markowitz,w_MaxSR)

windows()
barplot(t(wMV),beside=TRUE,ylab="Part (%)",ylim=c(0,0.4))
legend("topright",legend = c("PMVG","PT"), fill=c("black","gray"))

## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Portafolios Robustos

# Incertidumbre de Intervalo

# Calculo delta (distribución normal)
t <- nrow(retornos[,1])
alpha <- 0.05
znorm <- qnorm((1-alpha/2), 0, 1)
(delta <- znorm*sigma/sqrt(t))

# Calculo lambda: COef. Sharpe
rpSharpe <- mu%*%w_MaxSR
sigmapSharpe <- sqrt(t(w_MaxSR)%*%cov%*%w_MaxSR)
(lambda <- (rpSharpe-rf)/sigmapSharpe^2)

portfolioRobustBoxU <- function(mu, cov, delta){
    w <- Variable(length(mu))
    prob <- Problem(Maximize( t(w)%*%mu-t(abs(w))%*%delta-lambda*(quad_form(w, cov)) ),
                    constraints = list(w >= 0, ones(1,n)%*%w==1)) #sum(w) == 1
    result <- solve(prob)
    return(as.vector(result$getValue(w)))
}

w_RobustBox <- round(portfolioRobustBoxU(mu, cov, delta),4)
names(w_RobustBox) <- colnames(cov)
(w_RobustBox  <- cbind(w_RobustBox))

windows()
barplot(t(w_RobustBox), main="PRMVi",beside = TRUE,ylab="Part (%)",
        ylim=c(0,0.3))

wcomp <- cbind(w_RobustBox,w_MaxSR)

windows()
barplot(t(wcomp),beside=TRUE,ylab="Part (%)",ylim=c(0,0.4))
legend("topright",legend = c("PRMVi","PT"), fill=c("black","gray"))


## ----------------------------------------------------------------

# Ejemplo 13.2
# Incertidumbre Elipsoidal 

(delta2 = sqrt(qchisq(alpha,n)))

portfolioRobustEllipsoidU <- function(mu, cov, delta) {
    S12 <- chol(cov)  # t(S12) %*% S12 = Sigma
    w <- Variable(length(mu))
    prob <- Problem(Maximize( t(w)%*%mu-lambda*(quad_form(w, cov))-delta2*(norm2(S12%*%w))),
                    constraints = list(w >= 0, ones(1,n)%*%w==1))
    result <- solve(prob)
    return(as.vector(result$getValue(w)))
}

w_RobustEllipsoid <- round(portfolioRobustEllipsoidU(mu, cov, delta2),4)
names(w_RobustEllipsoid) <- colnames(cov)
(w_RobustEllipsoid  <- cbind(w_RobustEllipsoid))

# COmparación de resultados
windows()
barplot(t(w_RobustEllipsoid), main="PRMVe",beside = TRUE,ylab="Part (%)",
        ylim=c(0,0.3))

wcomp2 <- cbind(w_RobustEllipsoid,w_MaxSR)

windows()
barplot(t(wcomp2),beside=TRUE,ylab="Part (%)",ylim=c(0,0.4))
legend("topright",legend = c("PRMVe","PT"), fill=c("black","gray"))

## ----------------------------------------------------------------
## ----------------------------------------------------------------

# Evaluacion de desempeno 

# Dentro de muestra

# Retornos
(rp.PT <- t(w_MaxSR)%*%mu)
(rp.PRMVi <- t(w_RobustBox)%*%mu)
(rp.PRMVe <- t(w_RobustEllipsoid)%*%mu)
(mu.indice)

# Riesgo
(sigma.PT <- sqrt(t(w_MaxSR)%*%cov%*%w_MaxSR))
(sigma.PRMVi <- sqrt(t(w_RobustBox)%*%cov%*%w_RobustBox))
(sigma.PRMVe <- sqrt(t(w_RobustEllipsoid)%*%cov%*%w_RobustEllipsoid))
(sigma.indice)

# Coef. Sharpe
(S.PT <- (rp.PT-rf)/sigma.PT)
(S.PRMVi <- (rp.PRMVi-rf)/sigma.PRMVi)
(S.PRMVe <- (rp.PRMVe-rf)/sigma.PRMVe)
(S.Indice <- (mu.indice-rf)/sigma.indice)

# Desempeño historico - Valor del portafolio
valor <- 100 
t <- nrow(retornos)
n <- ncol(retornos)

# Retornos historicos

# Portafolio PT
rppt <- retornos%*%w_MaxSR
# PRMVi
rprmvi <- retornos%*%w_RobustBox
# PRMVe
rprmve <- retornos%*%w_RobustEllipsoid

# Valor del portafolio PT
port.pt <- matrix(0, nrow=t)
port.pt[1] <- valor
for(i in 2:t){
    port.pt[i] <- port.pt[i-1]*exp(rppt[i-1])
}

# Valor del PRMVi
port.rmvi <- matrix(0, nrow=t)
port.rmvi[1] <- valor
for(i in 2:t){
    port.rmvi[i] <- port.rmvi[i-1]*exp(rprmvi[i-1])
}

# Valor del PRMVe
port.rmve <- matrix(0, nrow=t)
port.rmve[1] <- valor
for(i in 2:t){
    port.rmve[i] <- port.rmve[i-1]*exp(rprmve[i-1])
}

# Valor del benchmark
v.benchmark <- matrix(0, nrow=t)
v.benchmark[1] <- valor
for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.indice[i-1])
}

Performance <- cbind(port.rmvi,port.rmve,port.pt,v.benchmark)
Performance <- ts(Performance,start=2016, frequency=12)
colnames(Performance) <- c("PRMVi","PRMVe","PT","Benchmark")
                     
windows()
plot(Performance[,1], type='l',main="Evaluación de desempeño",ylim=c(100,400), 
     ylab="Valor del portafolio", xlab="Tiempo (años)",bty="L")
lines(Performance[,2], col='darkblue')
lines(Performance[,3],lty=2)
lines(Performance[,4], col='darkgray')
legend("topleft",c("PRMVi","PRMVe","PT","S&P 500"),
       lty =c(1,1,2,1), col=c("black","darkblue","black","darkgray"))

## ----------------------------------------------------------------

# Desempeño Out-sample

fechaifm <- '2020-12-01'
fechaffm <- '2021-12-31'

preciosfm <- precios(activos,fechaifm,fechaffm,periodicidad)
retornosfm <- diff(log(preciosfm))[-1,]
indicefm <- precios(indice,fechaifm,fechaffm,periodicidad)
r.indicefm <- diff(log(indicefm))[-1,]

valor <- 100 
t <- nrow(retornosfm)+1
n <- ncol(retornosfm)

## Retornos historicos
rpptfm <- retornosfm%*%w_MaxSR
rprmvifm <- retornosfm%*%w_RobustBox
rprmvefm <- retornosfm%*%w_RobustEllipsoid

## Retornos esperados
(mu.PTfm <- mean(rpptfm))
(mu.PRMVifm <- mean(rprmvifm))
(mu.PRMVefm <- mean(rprmvefm))
(mu.indicefm <- mean(r.indicefm))

# Riesgo
(sigma.PTfm <- sd(rpptfm))
(sigma.PRMVifm <- sd(rprmvifm))
(sigma.PRMVefm <- sd(rprmvefm))
(sigma.indicefm <- sd(r.indicefm))

# Coef. Sharpe
(S.PTfm <- (mu.PTfm-rf)/sigma.PTfm)
(S.PRMVifm <- (mu.PRMVifm-rf)/sigma.PRMVifm)
(S.PRMVefm <- (mu.PRMVefm-rf)/sigma.PRMVefm)
(S.Indicefm <- (mu.indicefm-rf)/sigma.indicefm)

# Valor del portafolio PT
port.ptfm <- matrix(0, nrow=t)
port.ptfm[1] <- valor
for(i in 2:t){
    port.ptfm[i] <- port.ptfm[i-1]*exp(rpptfm[i-1])
}

# Valor del PRMVi
port.rmvifm <- matrix(0, nrow=t)
port.rmvifm[1] <- valor
for(i in 2:t){
    port.rmvifm[i] <- port.rmvifm[i-1]*exp(rprmvifm[i-1])
}

# Valor del PRMVe
port.rmvefm <- matrix(0, nrow=t)
port.rmvefm[1] <- valor
for(i in 2:t){
    port.rmvefm[i] <- port.rmvefm[i-1]*exp(rprmvefm[i-1])
}

# Benchmark
v.benchmarkfm <- matrix(0, nrow=t)
v.benchmarkfm[1] <- valor
for(i in 2:t){
    v.benchmarkfm[i] <- v.benchmarkfm[i-1]*exp(r.indicefm[i-1])
}

Performancefm <- cbind(port.rmvifm,port.rmvefm,port.ptfm,v.benchmarkfm)
Performancefm <- ts(Performancefm,start=c(2021,1), frequency=12)
colnames(Performancefm) <- c("PRMVi","PRMVe","PT","Benchmark")

windows()
plot(Performancefm[,1], type='l',main="Evaluación de desempeño", 
     ylab="Valor del portafolio", xlab="Tiempo (meses)",bty="L") #,ylim=c(100,400)
lines(Performancefm[,2], col='darkblue')
lines(Performancefm[,3],lty=2)
lines(Performancefm[,4], col='darkgray')
legend("topleft",c("PRMVi","PRMVe","PT","S&P 500"),
       lty =c(1,1,2,1), col=c("black","darkblue","black","darkgray"))


#----------------------------------------------------------------------
#----------------------------------------------------------------------

# Comparación de resultados
# Ventana movil

Sharpe_ratio <- function(w, mu, cov) {
    return(t(w) %*% mu / sqrt(t(w) %*% cov %*% w))
    }

Sharpe_Markowitz <- NULL
Sharpe_RobustEllip <- NULL

nport <- 100
w_all_Markowitz <- matrix(0,ncol=n,nrow=nport)
w_all_Markowitz_RobEllip <- matrix(0,ncol=n,nrow=nport)

set.seed(123)
for (i in 1:nport) {
    ret_rand <- mvrnorm(t, mu, cov)
    mu_rand <- colMeans(ret_rand)
    cov_rand <- cov(ret_rand)
    
    w_Markowitz <- cbind(portfolioMarkowitz(mu_rand,cov_rand))
    w_RMVe <- cbind(portfolioRobustEllipsoidU(mu_rand, cov_rand, delta2))
    
    w_all_Markowitz[i,] <- rbind(w_Markowitz)
    w_all_Markowitz_RobEllip[i,] <- rbind(w_RMVe)
    
    Sharpe_Markowitz[i] <- apply(w_Markowitz,MARGIN=2, FUN=Sharpe_ratio,mu_rand,cov_rand)
    Sharpe_RobustEllip[i] <- apply(w_RMVe,MARGIN=2,FUN=Sharpe_ratio,mu_rand,cov_rand)
}

(meanMV <- mean(Sharpe_Markowitz))
(meanRMVe <- mean(Sharpe_RobustEllip))
(SR_fix <- Sharpe_ratio(w_Markowitz, mu, cov))

windows()
plot(Sharpe_rand, col="darkgreen",pch=19)
points(Sharpe_RobustEllip, col="black",pch=19)
abline(h = SR_fix)
legend("topleft",c("MV","RobustEllip"), 
       col=c("darkgreen","black"),pch=c(19,19)) 
text((nport-2),SR_fix,labels = paste("Sharpe = ",round(SR_fix,4)),pos = 3)


#----------------------------------------------------------------------
#----------------------------------------------------------------------
    











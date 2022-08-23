
#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestion Moderna de portafolio
# Capitulo 15: Maxima diversificacion
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

(mu <- colMeans(retornos))
(cov <- cov(retornos))
var <- diag(cov)
(sigma <- sqrt(var))
n <- length(mu)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

# Ejemplo 15.1
## Portafolio Optimo de Markowitz: PMVG
unos <- rep(1,n)
cov_inv_1 <- solve(cov, unos) 
(wpmvg <- round((1/as.numeric(unos %*% cov_inv_1)) * cov_inv_1,4))
(rpmvg <- round(mu%*%wpmvg,4))
(sigmapmvg <- round(sqrt(t(wpmvg)%*%cov%*%wpmvg),4))
(drpmvg <- round(sigma%*%wpmvg/sigmapmvg,2))

## Portafolio Max DR

# Formulacion general

wpmd <- solve(cov, sigma)/as.numeric(t(sigma)%*%solve(cov)%*%sigma)
(wpmd = round(wpmd / sum(wpmd),4))
(rpmd <- round(mu%*%wpmd,4))
(sigmapmd <- round(sqrt(t(wpmd)%*%cov%*%wpmd),4))
(drpmd <- round(sigma%*%wpmd/sigmapmd,2))

# ------------------------------------------------------------------

# Portafolo Max DR sin cortos
# Uso de quadprog

lb = rep(0, n)
ub = rep(1, n)

Alb = -lb %*% matrix(1, 1, n) + diag(n) # Equiv: diag(n)
Aub = ub %*% matrix(1, 1, n) - diag(n)

Amat = t(rbind(sigma, Alb, Aub))
bvec = c(1,rep(0,2*n))
wpmd2 = solve.QP(cov,c(rep(0,n)),Amat,bvec,1)

wpmd2 = wpmd2$solution
(wpmd2 = wpmd2 / sum(wpmd2))
names(wpmd2) <- activos
sum(wpmd2)

# ------------------------------------------------------------

# Comparacion de resultados con el PMVG
wcomp <- rbind(wpmd,wpmvg)

windows()
barplot(wcomp,beside=TRUE,ylab="Part (%)",ylim=c(0,0.25))
legend("topright",legend = c("MDR","PMVG"), fill=c("black","gray"))

#---------------------------------------------------------------------

# Formulacion MDR Naive
sigmadiag <- diag(sigma^2)
wpmdn <- solve(sigmadiag, sigma)/as.numeric(t(sigma)%*%solve(sigmadiag)%*%sigma)
(wpmdn = round(wpmdn / sum(wpmdn),4))
(rpmdn <- round(mu%*%wpmdn,4))
(sigmapmdn <- round(sqrt(t(wpmdn)%*%cov%*%wpmdn),4))
(drpmdn <- round(sigma%*%wpmdn/sigmapmdn,2))

# Portafolio NRP
wrpnaive <- 1/sigma/(sum(1/sigma))

wcomp2 <- rbind(wpmdn,wrpnaive)

windows()
barplot(wcomp2,beside=TRUE,ylab="Part (%)",ylim=c(0,0.2))
legend("topright",legend = c("MDR","NRP"), fill=c("black","gray"))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

# Evaluacion de desempeno 

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
rpmvhist <- retornos%*%wpmvg
(mu.rpmvg <- mean(rpmvhist))
(sigma.pmvg <- sd(rpmvhist))
(sharpe.pmvg <- mu.rpmvg/sigma.pmvg)

# Portafolio DR
rpmdhist <- retornos%*%wpmd
(mu.rpmd <- mean(rpmdhist))
(sigma.md <- sd(rpmdhist))
(sharpe.md <- mu.rpmd/sigma.md)

# Valor del portafolio
# PMVG
port.mv <- matrix(0, nrow=t)
port.mv[1] <- valor
for(i in 2:t){
    port.mv[i] <- port.mv[i-1]*exp(rpmv[i-1])
}

# Paridad de Riesgo
port.md <- matrix(0, nrow=t)
port.md[1] <- valor
for(i in 2:t){
    port.md[i] <- port.md[i-1]*exp(rpmdhist[i-1])
}

# Benchmark
v.benchmark <- matrix(0, nrow=t)
v.benchmark[1] <- valor
for(i in 2:t){
    v.benchmark[i] <- v.benchmark[i-1]*exp(r.indice[i-1])
}

Performance <- cbind(port.md,port.mv,v.benchmark)
Performance <- ts(Performance,start=2016, frequency=12)
colnames(Performance) <- c("MDR","PMVG","Benchmark")

windows()
plot(Performance[,1], type='l',main="Evaluacion de desempeño", 
     ylab="Valor del portafolio", xlab="Tiempo (años)",bty="L")
lines(Performance[,2], lty=2)
lines(Performance[,3], col='darkgray')
legend("topleft",c("MDR","PMVG","S&P 500"),
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

# MDR
rpmdfm <- retornosfm%*%wpmd
(mu.rpmdfm <- mean(rpmdfm))
(sigma.mdfm <- sd(rpmdfm))
(sharpe.mdfm <- mu.rpmdfm/sigma.mdfm)

# Valor del portafolio
# PMVG
port.mvfm <- matrix(0, nrow=t)
port.mvfm[1] <- valor
for(i in 2:t){
    port.mvfm[i] <- port.mvfm[i-1]*exp(rpmvfm[i-1])
}

# MDR
port.mdfm <- matrix(0, nrow=t)
port.mdfm[1] <- valor
for(i in 2:t){
    port.mdfm[i] <- port.mdfm[i-1]*exp(rpmdfm[i-1])
}

# Benchmark
v.benchmarkfm <- matrix(0, nrow=t)
v.benchmarkfm[1] <- valor
for(i in 2:t){
    v.benchmarkfm[i] <- v.benchmarkfm[i-1]*exp(r.indicefm[i-1])
}

Performancefm <- cbind(port.mdfm,port.mvfm,v.benchmarkfm)
Performancefm <- ts(Performancefm,start=c(2021,1), frequency=12)
colnames(Performancefm) <- c("MDR","PMVG","Benchmark")

windows()
plot(Performancefm[,1], type='l',main="Evaluacion de desempeno", bty="L",
     ylab="Valor del portafolio", xlab="Tiempo (meses)", ylim=c(min(Performancefm), max(Performancefm)))
lines(Performancefm[,2], lty=2)
lines(Performancefm[,3], col='darkgray')
legend("topleft",c("MDR","PMVG","S&P 500"),
       lty =c(1,2,1), col=c("black","black","darkgray"))





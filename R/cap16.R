
#-------------------------------------------------------------
#-------------------------------------------------------------
# Gestion Moderna de portafolio
# Capitulo 16: Portafolios sostenibles y criterios ASG
# Copyright 2022
#-------------------------------------------------------------
#-------------------------------------------------------------
#
# Notas:
# 1. Requiere cargar primero la funcion para importar precios  
#    de Yahoo Finance: "func_precios"
#-------------------------------------------------------------
#-------------------------------------------------------------

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

# ------------------------------------------------------------
## Scraping indicadores ASG de Yahoo Finance 
## Fuente: Sustainalytics - MorningStar
## -----------------------------------------------------------
#
# Fuente: https://christophscheuch.github.io/post/r-stuff/scraping-esg-data/
## -----------------------------------------------------------

library(dplyr)
library(tidyverse) 
library(Rcpp)  

asg_func <- function(symbol){
    scrap_date <- fechaf #Sys.time()
    url <- str_c("https://finance.yahoo.com/quote/",symbol, "/sustainability?p=", symbol)
    page <- GET(url) %>% read_html()
    asg_score <- page %>% 
        html_node(xpath = '//*[@id="Col1-0-Sustainability-Proxy"]/section/div[1]/div/div[1]/div/div[2]/div[1]') %>%
        html_text() %>% parse_number()
    return(asg_score)
}

iasg <- lapply(activos,asg_func)
(iasg <- as.numeric(asg_data))

(datos_resumen <- round(rbind(mu,sigma,iasg),4))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

## FE de Markowitz
unos <- rep(1,n)
x <- t(mu)%*%solve(cov)%*%mu
y <- t(mu)%*%solve(cov)%*%unos
z <- t(unos)%*%solve(cov)%*%unos
d <- x*z - y*y
g <- (solve(cov,unos)%*%x-solve(cov,mu)%*%y)%*%solve(d)
h <- (solve(cov,mu)%*%z-solve(cov,unos)%*%y)%*%solve(d)

nport <- 100
j <- seq(min(mu),max(mu), length=nport) 
wpo <- matrix(c(0), ncol=n, nrow=nport) 
rpo <- matrix(c(0), nrow=nport)
sigmapo <- matrix(c(0), nrow=nport)
for(i in 1:nport){
    wj <- g + h*j[i] 
    wpo[i,] <- t(wj)
    rpo[i,] <- t(wj)%*%mu
    sigmapo[i,] <- sqrt(t(wj)%*%cov%*%wj)
}

## PMVG
cov_inv_1 <- solve(cov, ones) 
(wpmvg <- round((1/as.numeric(ones %*% cov_inv_1)) * cov_inv_1,4))
(rpmvg <- round(mu%*%wpmvg,4))
(sigmapmvg <- round(sqrt(t(wpmvg)%*%cov%*%wpmvg),4))
(ipmvg <- t(wpmvg)%*%iasg)

# ------------------------------------------------------------------

# Modelo MV-ESG

unos <- rep(1,n)
x <- t(mu)%*%solve(cov)%*%mu
y <- t(mu)%*%solve(cov)%*%unos
z <- t(unos)%*%solve(cov)%*%unos

# Extension ESG
p <- t(mu)%*%solve(cov)%*%iasg
q <- t(unos)%*%solve(cov)%*%iasg
r <- t(iasg)%*%solve(cov)%*%iasg

# Calculo d ajustado
d <- (x*z*r+y*p*q+y*p*q-z*p^2-x*q^2-r*y^2)

# Vectores
g <- (solve(cov,mu)%*%(p*q)-solve(cov,mu)%*%(y*r)+solve(cov,unos)%*%(x*r)-solve(cov,unos)%*%(p^2)+
          solve(cov,iasg)%*%(p*y)-solve(cov,iasg)%*%(x*q))%*%(1/d)
h <- (solve(cov,mu)%*%(r*z)-solve(cov,mu)%*%(q^2)+solve(cov,unos)%*%(p*q)-solve(cov,unos)%*%(y*r)+
          solve(cov,iasg)%*%(y*q)-solve(cov,iasg)%*%(p*z))%*%(1/d)
f <- (solve(cov,mu)%*%(y*q)-solve(cov,mu)%*%(z*p)+solve(cov,unos)%*%(y*p)-solve(cov,unos)%*%(q*x)+
          solve(cov,iasg)%*%(z*x)-solve(cov,iasg)%*%(y^2))%*%(1/d)

# -----------------------------------------------------------------
# Pesos optimos casos especiales

# Caso: mean(ASG_score), mean(Ri)

(ip <- mean(iasg))
(rp <- mean(mu))
(wpmv1 <- g+h*rp+f*ip)

pesos <- cbind(wpmv1,wpmvg)

windows()
barplot(t(pesos),beside=TRUE,ylab="Part (%)",ylim=c(0,0.25))
legend("topright",legend = c("MV-ASG","PMVG"), fill=c("black","gray"))

# Caso: Ip=Ip(MVG), Rp=Rp(MVG)

(ipMVG <- c(t(wpmvg)%*%iasg))
(rpMVG <- c(rpmvg))
(wpmv2 <- g+h*rpMVG+f*ipMVG)

pesos2 <- cbind(wpmv2,wpmvg)

windows()
barplot(t(pesos2),beside=TRUE,ylab="Part (%)",ylim=c(0,0.25))
legend("topright",legend = c("MV-ASG","PMVG"), fill=c("black","gray"))

#--------------------------------------------------------------------
#--------------------------------------------------------------------

# Construccion FE-ASG - Para: iasg < ipmvg

nport <- 100
rpj <- seq(min(mu),max(mu),length=nport)
ipj <- seq(min(iasg),ipMVG,length=5)

wpoASG <- matrix(0,ncol=n,nrow=nport)
sigmapoASG <- matrix(0,ncol=nport,nrow=nport)
rpoASG <- matrix(0,ncol=nport,nrow=nport)
ipoASG <- matrix(0,ncol=5,nrow=nport) 

for(i in 1:nport){
    for(j in 1:5){
        (w <- g+h*rpj[i]+f*ipj[j])
        wpoASG[i,] <- t(w)
        sigmapoASG[i,j] <- sqrt(t(w)%*%cov%*%w)
        rpoASG[i,j] <- t(w)%*%mu
        ipoASG[i,j] <- t(w)%*%iasg
    }
}

windows()
plot(sigma,mu,main="Plano Riesgo-Retorno",xlim=c(0,max(sigma)*1.2),
     ylim=c(0.01,max(mu)),xlab="Riesgo",ylab="Retorno esperado",
     pch=19,bty="L")
points(sigmapoASG[,1],rpoASG[,1],type="l",lty=6)
points(sigmapoASG[,2],rpoASG[,2],type="l",lty=2)
points(sigmapoASG[,3],rpoASG[,3],type="l",lty=3)
points(sigmapoASG[,4],rpoASG[,4],type="l",lty=4)
points(sigmapoASG[,5],rpoASG[,5],type="l",lty=5)
lines(sigmapo,rpo,lwd=3)
text(sigma,mu,labels = activos,cex=0.6,pos=4)
points(sigmapmvg,rpmvg,pch=19)
text(sigmapmvg,rpmvg,labels = "PMVG",cex=0.6,pos=2)
legend("topleft",legend=c("FE","FE-ASG"),lty = c(1,2), lwd=c(2,1),
       col=c("black","black"))

#------------------------------------------------------------------

# Construccion FE-ASG - Para: iasg > ipmvg

nport <- 100
rpj <- seq(min(mu),max(mu),length=nport)
ipj <- seq(ipMVG,max(iasg),length=5)

wpoASG <- matrix(0,ncol=n,nrow=nport)
sigmapoASG <- matrix(0,ncol=nport,nrow=nport)
rpoASG <- matrix(0,ncol=nport,nrow=nport)
ipoASG <- matrix(0,ncol=5,nrow=nport) 

for(i in 1:nport){
    for(j in 1:5){
        (w <- g+h*rpj[i]+f*ipj[j])
        wpoASG[i,] <- t(w)
        sigmapoASG[i,j] <- sqrt(t(w)%*%cov%*%w)
        rpoASG[i,j] <- t(w)%*%mu
        ipoASG[i,j] <- t(w)%*%iasg
    }
}

# Graficos
windows()
plot(sigma,mu,main="Plano Riesgo-Retorno",xlim=c(0,max(sigma)*1.2),
     ylim=c(0.01,max(mu)), xlab="Riesgo",ylab="Retorno esperado",
     pch=19,bty="L")
points(sigmapoASG[,1],rpoASG[,1],type="l",lty=6)
points(sigmapoASG[,2],rpoASG[,2],type="l",lty=2)
points(sigmapoASG[,3],rpoASG[,3],type="l",lty=3)
points(sigmapoASG[,4],rpoASG[,4],type="l",lty=4)
points(sigmapoASG[,5],rpoASG[,5],type="l",lty=5)
lines(sigmapo,rpo,lwd=3)
text(sigma,mu,labels = activos,cex=0.6,pos=4)
points(sigmapmvg,rpmvg,pch=19)
text(sigmapmvg,rpmvg,labels = "PMVG",cex=0.6,pos=2)
legend("topleft",legend=c("FE","FE-ASG"),lty = c(1,2), lwd=c(2,1),
       col=c("black","black"))

# ----------------------------------------------------------------
#------------------------------------------------------------------

# Construccion Superficie SE-ASG 

nport <- 100
rpj <- seq(min(mu),max(mu),length=nport)
ipj <- seq(min(iasg),max(iasg),length=nport)

wpoASG <- matrix(0,ncol=n,nrow=nport)
sigmapoASG <- matrix(0,ncol=nport,nrow=nport)
rpoASG <- matrix(0,ncol=nport,nrow=nport)
ipoASG <- matrix(0,ncol=nport,nrow=nport) 

for(i in 1:nport){
    for(j in 1:nport){
        (w <- g+h*rpj[i]+f*ipj[j])
        wpoASG[i,] <- t(w)
        sigmapoASG[i,j] <- sqrt(t(w)%*%cov%*%w)
        rpoASG[i,j] <- t(w)%*%mu
        ipoASG[i,j] <- t(w)%*%iasg
    }
}

# Graficos
windows()
plot(sigma,mu,main="Plano Riesgo-Retorno",xlim=c(0,max(sigma)*1.2),
     ylim=c(0.01,max(mu)), xlab="Riesgo",ylab="Retorno esperado",
     pch=19,bty="L")
points(sigmapoASG,rpoASG,col="gray")
lines(sigmapo,rpo,lwd=3)
text(sigma,mu,labels = activos,cex=0.6,pos=4)
points(sigmapmvg,rpmvg,pch=19)
text(sigmapmvg,rpmvg,labels = "PMVG",cex=0.6,pos=2)
legend("topleft",legend=c("FE","FE-ASG"),lty = c(1,1), lwd=c(2,1),
       col=c("black","gray"))


windows()
plot(sigma,iasg,main="Plano Riesgo-Retorno",xlim=c(0,max(sigma)*1.2),ylim=c(0,max(iasg)),
     xlab="Riesgo",ylab="ESG score",pch=19)
points(sigmapoASG,ipoASG,col="gray")
text(sigma,iasg,labels = activos,cex=0.6,pos=4)
points(sigmapmvg,ipmvg,pch=19)
text(sigmapmvg,ipmvg,labels = "PMVG",cex=0.7,pos=2)


# Grafico 3D
library(plotly)
# Plot
p <- plot_ly(x=sigmapoASG,y=rpoASG,z =ipoASG, type = "surface")
p %>% add_lines(x=sigmapo,y=rpo)
p %>% layout(scene = list(
    xaxis = list(title = "Riesgo"),
    yaxis = list(title = "Retorno esperado"),
    zaxis = list(title = "ESG score")
)
)


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





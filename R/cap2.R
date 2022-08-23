
#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestion Moderna de portafolio
# Capitulo 2: Modelo MV
# Copyright 2022
#-----------------------------------------------------------
#-----------------------------------------------------------
#
# Notas:
# 1. Requiere cargar primero la funcion para importar precios  
#    de Yahoo Finance: "func_precios"
#-----------------------------------------------------------
#-----------------------------------------------------------

## Ejemplos 2.1 y 2.2

rm(list=ls())
library(zoo)
library(xts)
library(quantmod)
library(quadprog)

## --------------------------------------------------------
## Caso de 2 activos: ejemplo 2.1
## --------------------------------------------------------

# Se importa la informacion de precios para las acciones 

# Activos: AAPL y AMZN
fechai <- "2009-12-01"
fechaf <- "2020-12-31"
periodicidad <- "monthly" 
activos <- c("AAPL","AMZN")
precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,] 

(mu <- colMeans(retornos)*12)
cov <- cov(retornos)*12
var <- diag(cov)
(sigma <- sqrt(var))
(rho <- cor(retornos$AAPL,retornos$AMZN))
n <- length(mu)

# Calculo de los pesos optimos

(w1 <- (sigma[2]^2 - rho*sigma[1]*sigma[2])/(sigma[1]^2+sigma[2]^2 - 2*(rho*sigma[1]*sigma[2])))
(w2 <- (sigma[1]^2 - rho*sigma[1]*sigma[2])/(sigma[1]^2+sigma[2]^2 - 2*(rho*sigma[1]*sigma[2])))

# Portafolio  optimo
wp <- c(w1,w2)
(rp <- t(wp)%*%mu)
(sigmap <- sqrt(t(wp)%*%cov%*%wp))

windows()
par(mfrow=c(1,1),bty="L")
plot(sigma,mu,main="Plano Riesgo-Retorno",xlim=c(0.2,max(sigma)*1.1),
     ylim=c(min(mu)*0.9,max(mu)*1.1),
     xlab="Riesgo",ylab="Retorno esperado")
points(sigmap,rp,pch=19)
text(sigma,mu,labels=activos,cex=0.8,pos=4)
text(sigmap,rp,labels="P",pos=2)
abline(v=min(sigma),lty=2,col="gray")
arrows(min(sigma)*0.99, rp, sigmap*1.04, rp)
text(sigmap*1.1,rp*0.99,labels="Diversificación",pos=1)

## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Cambios en la correlacion de los activos
## Ejemplo 2.2

(rho1 <- 0.7 )
(w11 <- (sigma[2]^2 - rho1*sigma[1]*sigma[2])/(sigma[1]^2+sigma[2]^2 - 2*(rho1*sigma[1]*sigma[2])))
(w21 <- (sigma[1]^2 - rho1*sigma[1]*sigma[2])/(sigma[1]^2+sigma[2]^2 - 2*(rho1*sigma[1]*sigma[2])))
(sigmap1 <- sqrt(sigma[1]^2*w11^2+sigma[2]^2*w21^2+2*(rho1*w11*w21*sigma[1]*sigma[2])))

(rho2 <- -0.3)
(w12 <- (sigma[2]^2 - rho2*sigma[1]*sigma[2])/(sigma[1]^2+sigma[2]^2 - 2*(rho2*sigma[1]*sigma[2])))
(w22 <- (sigma[1]^2 - rho2*sigma[1]*sigma[2])/(sigma[1]^2+sigma[2]^2 - 2*(rho2*sigma[1]*sigma[2])))
(sigmap2 <- sqrt(sigma[1]^2*w12^2+sigma[2]^2*w22^2+2*(rho2*w12*w22*sigma[1]*sigma[2])))

windows()
par(mfrow=c(1,1),bty="L")
plot(sigma,mu,main=" ",xlim=c(0.12,max(sigma)*1.1),
     ylim=c(min(mu)*0.9,max(mu)*1.1),
     xlab="Riesgo",ylab="Retorno esperado")
points(sigmap,rp,pch=19)
text(sigma,mu,labels=activos,cex=0.8,pos=4)
text(sigmap*0.99,rp,labels=expression(rho),pos=3,cex=0.7)
text(sigmap*1.035,rp,labels=paste("=",round(rho,3)),pos=3,cex=0.7)
arrows(min(sigma)*0.99, mu[2], sigmap, mu[2],angle=10)
points(sigmap1,rp,pch=19)
text(sigmap1,rp,labels=expression(rho),pos=3,cex=0.7)
text(sigmap1*1.03,rp,labels=paste("=",rho1),pos=3,cex=0.7)
arrows(min(sigma)*0.99, mu[2]*1.01, sigmap1, mu[2]*1.01,angle=10)
points(sigmap2,rp,pch=19)
text(sigmap2,rp,labels=expression(rho),pos=3,cex=0.7)
text(sigmap2*1.042,rp,labels=paste("=",rho2),pos=3,cex=0.7)
arrows(min(sigma)*0.99, mu[2]*0.99, sigmap2, mu[2]*0.99,angle=10)
abline(v=min(sigma),lty=2,col="gray")

## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Caso de n activos
## Ejemplos 2.3, 2.4 y 2.5

## Acciones seleccionadas
fechai <- "2009-12-01"
fechaf <- "2021-12-31"
periodicidad <- "monthly" 
activos <- c("AAPL","AMZN","GOOGL","MSFT")

precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]

(mu <- colMeans(retornos)*12)
(cov <- cov(retornos)*12)
var <- diag(cov)
(sigma <- sqrt(var))
(corr <- cor(retornos))

## Portafolio Optimo de Markowitz
n <- ncol(retornos) 
ones <- rep(1,n)
x <- t(mu)%*%solve(cov)%*%mu
y <- t(mu)%*%solve(cov)%*%ones
z <- t(ones)%*%solve(cov)%*%ones
d <- x*z - y*y
g <- (solve(cov,ones)%*%x-solve(cov,mu)%*%y)%*%solve(d)
h <- (solve(cov,mu)%*%z-solve(cov,ones)%*%y)%*%solve(d)

# Un solo portafolio con retorno: Rp=0.25
rpobj <- 0.25
(wpobj <- g + h*rpobj)
rownames(wpobj) <- activos
(sigmapobj <- sqrt(t(wpobj)%*%cov%*%wpobj))

# Vector de pesos
windows()
barplot(t(wpobj))

# Para toda la Frontera Efieciente
rpmin <- min(mu)
rpmax <- max(mu)
nport <- 100

j <- seq(rpmin,rpmax, length=nport) 
wpo <- matrix(c(0), ncol=n, nrow=nport) 
rpo <- matrix(c(0), nrow=nport)
sigmapo <- matrix(c(0), nrow=nport)
wj <- 0

for(i in 1:nport){
    wj <- g + h*j[i] 
    wpo[i,] <- t(wj)
    rpo[i,] <- t(wj)%*%mu
    sigmapo[i,] <- sqrt(t(wj)%*%cov%*%wj)
}

# PMVG
cov_inv_1 <- solve(cov, ones) 
(wpmvg <- (1/as.numeric(ones %*% cov_inv_1)) * cov_inv_1)
(rpmvg <- mu%*%wpmvg)
(sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg))

# Frontera eficiente: FE
# Para la Frontera Efieciente con pendiente positiva
rpmin <- rpmvg
rpmax <- max(mu)
nport <- 100

j <- seq(rpmin,rpmax, length=nport) 
wpo2 <- matrix(c(0), ncol=n, nrow=nport) 
rpo2 <- matrix(c(0), nrow=nport)
sigmapo2 <- matrix(c(0), nrow=nport)
wj <- 0

for(i in 1:nport){
    wj <- g + h*j[i] 
    wpo2[i,] <- t(wj)
    rpo2[i,] <- t(wj)%*%mu
    sigmapo2[i,] <- sqrt(t(wj)%*%cov%*%wj)
}

windows()
plot(sigma,mu, xlim=c(0.15,max(sigma)*1.1), ylim=c(0.16,max(mu)),cex=0.6,
     xlab="Riesgo",ylab="Retorno esperado",font = 6,bty="L")
points(sigmapo,rpo,type='l',lty=2)
points(sigmapo2,rpo2,type='l',lwd=2)
points(sigmapmvg,rpmvg,pch=19)
text(sigma,mu,labels=activos,pos = 4, cex=0.7,font = 3)
text(sigmapmvg,rpmvg,labels="Pmvg",pos = 2, cex=0.7,font = 3)

windows()
par(mfrow=c(1,4))
barplot(t(wpo[,1]),ylim=c(-0.1,0.6), main="AAPL") #
barplot(t(wpo[,2]),ylim=c(-0.1,0.5), main="AMZN")
barplot(t(wpo[,3]),ylim=c(-0.2,0.8), main="GOOG")
barplot(t(wpo[,4]),ylim=c(0,0.6), main="MSFT")

# ------------------------------------------------------
# ------------------------------------------------------

# Ejemplo 2.6
# Solucion del portafolio tangente: maximo Sharpe

# Formulacion 1
rf <- 0.0
r <- mu -rf # Excesos de retorno

# SOlucion: Z = V^-1 x Er
Z <- solve(cov,r)
sumaZ <- sum(Z)
(wpt <- Z/sumaZ)
round(wpt,4)
(rpt <- t(wpt)%*%mu)
(sigmapt <- sqrt(t(wpt)%*%cov%*%wpt))

# Construccion LMC
wpc <- seq(0,1.5,length=100)
rpc <- matrix(0,nrow=length(wpc))
sigmapc <- matrix(0,nrow=length(wpc))

for(i in 1:length(wpc)){
    rpc[i] <- wpc[i]*rpt+(1-wpc[i])*rf
    sigmapc[i]<-wpc[i]*sigmapt
}

# Plano Riesgo-Retorno
windows()
plot(sigma,mu, xlim=c(0.15,max(sigma)*1.1), ylim=c(0.15,max(mu)),
     cex=0.6,xlab="Riesgo",ylab="Retorno esperado",bty="L")
points(sigmapo,rpo,type='l')
points(sigmapmvg,rpmvg,pch=19,cex=1)
points(sigmapt,rpt,pch=19,cex=1)
text(sigma,mu,labels=activos,pos = 4, cex=0.7)
text(sigmapmvg,rpmvg,labels="PMVG",pos = 2, cex=0.7)
text(sigmapt,rpt,labels="T",pos = 2, cex=1)
points(sigmapc,rpc,col="black",lwd=2,type="l",lty=3)
legend("bottomright",legend = c("Frontera eficiente","LMC"),
       lty=c(1,2),cex=0.8)


# Figuras
# Precios hsitoricos
p1 <- ts(precios.hist, frequency = 12, start = c("2016"))
windows()
par(mfrow = c(2, 2))
plot(p1[,1],xlab="",ylab="AAPL")
plot(p1[,2],xlab="",ylab="AMZN")
plot(p1[,3],xlab="",ylab="MSFT")
plot(p1[,4],xlab="",ylab="GOOG")
logoing_func(r_logo, x=0.90, y=0.10, size=0.15)


# Pesos optimos del PT: maximo Sharpe
windows()
barplot(t(wpt), beside=TRUE,ylim=c(-0.6,1.20),cex.names = 0.75,
        ylab="Part. (%)")

#-------------------------------------------------------
#-------------------------------------------------------
# Forma general PT

(Er2 <- as.numeric(mu)%*%(1/y))
(wpt2 <- solve(cov)%*%Er2)
(rpt2 <- x/y)
(sigmapt2 <- sqrt(x/y^2))

#-------------------------------------------------------
#-------------------------------------------------------

# Ejemplo 2.7
# Solucion PMVG y PT sin cortos

# Programa QP - Restricciones en los cortos
library(quadprog)

# FE de Markowitz
nport <- 1000
Dmat <- cov*2
Amat <- cbind(mu,rep(1,n),diag(1,nrow=n))
dvec <- rep(0,n)

j <- seq(min(mu)*1.001,max(mu)*0.999,length=nport)

wpo3 <- matrix(0,nrow=nport,ncol=n)
sigmapo3 <- matrix(0,nrow=nport)
rpo3 <- matrix(0,nrow=nport)

for(i in 1:nport){
    bvec <- c(j[i],1,rep(0,n)) # restricciones de igualdad del sistema
    result <- solve.QP(Dmat, dvec, Amat, bvec, meq=2) #meq=restriccones de igualdad
    wj <- result[["solution"]]
    wpo3[i,] <- wj
    rpo3[i] <- mu%*%wj
    sigmapo3[i] <- sqrt(t(wj)%*%cov%*%wj)
}

#--------------------------------------------------------
# Portafolio Tangente 

sharpe_port <- (rpo3-rf)/sigmapo3
maxsharpe <- max(sharpe_port)
dif <- diff(sharpe_port)
length(dif[dif > 0])

windows()
plot(sharpe_port, type='l',xlab="Portafolio",ylab="Coef. Sharpe",
     ylim=c(min(sharpe_port),max(sharpe_port)*1.1),bty="L")
points(length(dif[dif > 0]),maxsharpe,pch=19)
text(length(dif[dif > 0]),maxsharpe,labels="Máx. Sharpe",pos = 3, cex=1)

tabla <- cbind(sharpe_port,rpo3,sigmapo3)
colnames(tabla) <- c("Sharpe","Retorno","Riesgo")
colnames(wpo) <- activos

tabla2 <- cbind(tabla,wpo3)

sort_tabla2 <- tabla2[order(-tabla2[,1]),]
port_maxsharpe <- round(sort_tabla2[1,],4)

wpt
wpt2 <- port_maxsharpe[4:length(port_maxsharpe)]
rpt2 <- port_maxsharpe[2]
sigmapt2 <- port_maxsharpe[3]

# Plano Riesgo-Retorno
windows()
plot(sigma,mu, xlim=c(min(sigma)*0.8,max(sigma)*1.1), 
     ylim=c(min(mu)*0.9,max(mu)*1.1),
     cex=0.6,xlab="Riesgo",ylab="Retorno esperado",bty="L")
points(sigmapo,rpo,type='l',lty=2,lwd=2)
points(sigmapt,rpt,pch=19)
points(sigmapt2,rpt2,pch=19)
text(sigma,mu,labels=activos,pos = 4, cex=0.7)
#text(sigmapt,rpt,labels="T",pos = 2, cex=1)
text(sigmapt2,rpt2,labels="T=T*",pos = 2, cex=1)
points(sigmapmvg,rpmvg,pch=19,cex=1)
text(sigmapmvg,rpmvg,labels="PMVG",pos = 2, cex=0.7)
points(sigmapo3,rpo3,col="black",lwd=2,type="l")
legend("bottomright",legend = c("FE con cortos","FE sin cortos"),
       lty=c(2,1),cex=0.9)

names(wpt2) <- activos 

windows()
barplot(t(wpt2), beside = TRUE, cex.names = 0.75,
        ylab="Part. (%)")

# ----------------------------------------------------------------
# ----------------------------------------------------------------

# Implementacion del algoritmo: Seccion 2.9

# ----------------------------------------------------------------

activos <- c("AAPL","ABT","AMZN","CAT","CSX","CSCO","GOOG","HD",
             "JNJ","MSFT","MCD","V") 

precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,] 

(mu <- colMeans(retornos)*12)
cov <- cov(retornos)*12
var <- diag(cov)
sigma <- sqrt(var)
corr <- cor(retornos)

## ---------------------------------------------------
## Portafolio Optimo de Markowitz (FE, PMVG y PT)
## Con cortos
## ---------------------------------------------------

n <- ncol(retornos) 
ones <- rep(1,n)
x <- t(mu)%*%solve(cov)%*%mu
y <- t(mu)%*%solve(cov)%*%ones
z <- t(ones)%*%solve(cov)%*%ones
d <- x*z - y*y
g <- (solve(cov,ones)%*%x-solve(cov,mu)%*%y)%*%solve(d)
h <- (solve(cov,mu)%*%z-solve(cov,ones)%*%y)%*%solve(d)

rpmin <- min(mu)
rpmax <- max(mu)
nport <- 100

j <- seq(rpmin,rpmax, length=nport) 
wpo <- matrix(c(0), ncol=n, nrow=nport) 
rpo <- matrix(c(0), nrow=nport)
sigmapo <- matrix(c(0), nrow=nport)
wj <- 0

for(i in 1:nport){
    wj <- g + h*j[i] 
    wpo[i,] <- t(wj)
    rpo[i,] <- t(wj)%*%mu
    sigmapo[i,] <- sqrt(t(wj)%*%cov%*%wj)
}

# PMVG
cov_inv_1 <- solve(cov, ones) 
wpmvg <- (1/as.numeric(ones %*% cov_inv_1)) * cov_inv_1
rpmvg <- mu%*%wpmvg
sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)

rf <- 0.0
r <- mu -rf 
Z <- solve(cov,r)
sumaZ <- sum(Z)
wpt <- Z/sumaZ
rpt <- t(wpt)%*%mu
sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)

windows()
plot(sigma,mu, xlim=c(min(sigma)*0.7,max(sigma)*1.1), 
     ylim=c(min(mu)*0.8,max(mu)*1.1),
     cex=0.6,xlab="Riesgo",ylab="Retorno esperado",bty="L")
points(sigmapo,rpo,type='l',lty=1,lwd=2)
points(sigmapt,rpt,pch=19)
text(sigma,mu,labels=activos,pos = 4, cex=0.7)
text(sigmapt,rpt,labels="T",pos = 2, cex=1)
points(sigmapmvg,rpmvg,pch=19,cex=1)
text(sigmapmvg,rpmvg,labels="PMVG",pos = 2, cex=0.7)

windows()
barplot(t(wpt), beside = TRUE, cex.names = 0.75,
        ylab="Part. (%)")

## --------------------------------------------------- 
## Sin cortos
## ---------------------------------------------------

library(quadprog)
nport <- 1000
Dmat <- cov*2
Amat <- cbind(mu,rep(1,n),diag(1,nrow=n))
dvec <- rep(0,n)

j <- seq(min(mu)*1.001,max(mu)*0.999,length=nport)

wpo2 <- matrix(0,nrow=nport,ncol=n)
sigmapo2 <- matrix(0,nrow=nport)
rpo2 <- matrix(0,nrow=nport)

for(i in 1:nport){
    bvec <- c(j[i],1,rep(0,n)) # restricciones de igualdad del sistema
    result <- solve.QP(Dmat, dvec, Amat, bvec, meq=2) #meq=restriccones de igualdad
    wj <- result[["solution"]]
    wpo2[i,] <- wj
    rpo2[i] <- mu%*%wj
    sigmapo2[i] <- sqrt(t(wj)%*%cov%*%wj)
}

sharpe_port <- (rpo2-rf)/sigmapo2
maxsharpe <- max(sharpe_port)
tabla <- cbind(sharpe_port,rpo2,sigmapo2)
colnames(wpo2) <- activos
tabla2 <- cbind(tabla,wpo2)
sort_tabla2 <- tabla2[order(-tabla2[,1]),]
port_maxsharpe <- round(sort_tabla2[1,],4)

wpt2 <- port_maxsharpe[4:length(port_maxsharpe)]
rpt2 <- port_maxsharpe[2]
sigmapt2 <- port_maxsharpe[3]

# Plano Riesgo-Retorno
windows()
plot(sigma,mu, xlim=c(min(sigma)*0.7,max(sigma)*1.2), ylim=c(min(mu)*0.8,max(mu)*1.2),
     cex=0.6,xlab="Riesgo",ylab="Retorno esperado",bty="L")
points(sigmapo,rpo,type='l',lty=2)
points(sigmapmvg,rpmvg,pch=19)
points(sigmapt,rpt,pch=19)
text(sigma,mu,labels=activos,pos = 4, cex=0.7)
text(sigmapmvg,rpmvg,labels="PMVG",pos = 2, cex=0.7)
text(sigmapt,rpt,labels="T",pos = 2, cex=1)
points(sigmapo2,rpo2,col="black",lwd=2,type="l")
points(sigmapt2,rpt2,pch=19)
text(sigmapt2,rpt2,labels="T*",pos = 4, cex=1)

windows()
barplot(t(wpt2),axisnames = TRUE, beside = TRUE, cex.names = 0.75,
        ylab="Part. (%)") #main="Sin cortos"



#--------------------------------------------------------
#------------------------------------
# Seccion 2.9
#------------------------------------
library(ggplot2)
library(PortfolioAnalytics)
library(corrplot)
library(ggcorrplot)

# COrrelaciones
windows()
chart.Correlation(retornos, histogram = F, cex=0.6)

corr2 <- round(corr, 1)

windows()
ggcorrplot(corr2, hc.order = TRUE, type = "lower",
           outline.col = "white",
           colors = c("gray", "darkgray","black")) +
    labs(title=" ", 
         subtitle=" ", 
         caption=" ", 
         color=NULL) +
    theme(legend.position="right")

data <- data.frame(mu,sigma)

p <- ggplot(data, aes(x = sigma, y = mu))  + geom_point(size = 3) +
    xlim(c(0, 0.1)) + 
    ylim(c(0, 0.04)) + 
    labs(y = "Retorno", 
         x = "Riesgo", 
         title = "") + 
    theme_classic()

windows()
p


#-------------------------------------------------------

# Ejemplo Limitaciones del modelo MV
# 
library(quadprog)
fechai <- "2009-12-01"
fechaf <- "2021-12-31"
periodicidad <- "monthly" 
activos <- c("AAPL","AMZN","GOOGL","MSFT")
precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]

n <- ncol(retornos)
t <- 60
tf <- nrow(retornos)
i <- c(t:tf)
nport <- 1000
wpt_roll <- matrix(0,nrow=length(i),ncol=n)
colnames(wpt_roll) <- activos

for(j in 1:length(i)){
    ret_roll <- retornos[1:i[j],] # retornos[i[j]-i[j]+1:i[j],]
    mu_roll <- colMeans(ret_roll)
    cov_roll <- cov(ret_roll)
    Dmat <- cov_roll*2
    Amat <- cbind(mu_roll,rep(1,n),diag(1,nrow=n))
    dvec <- rep(0,n)
    rp <- seq(min(mu_roll)*1.001,max(mu_roll)*0.999,length=nport)
    sigmapo_roll <- matrix(0,nrow=nport)
    rpo_roll <- matrix(0,nrow=nport)
    wpo_rollm <- matrix(0,nrow=nport, ncol=n)
    for(k in 1:nport){
        bvec <- c(rp[k],1,rep(0,n)) 
        result <- solve.QP(Dmat, dvec, Amat, bvec, meq=2) 
        wj <- result[["solution"]]
        rpo_roll[k] <- mu_roll%*%wj
        sigmapo_roll[k] <- sqrt(t(wj)%*%cov_roll%*%wj)
        wpo_rollm[k,] <- wj
    }
    sharpeport_roll <- (rpo_roll-rf)/sigmapo_roll
    tabla <- cbind(sharpeport_roll,rpo_roll,sigmapo_roll,wpo_rollm)
    sort_tabla <- tabla[order(-tabla[,1]),]
    port_maxsharpe <- round(sort_tabla[1,],4)
    wpt_roll[j,] <- port_maxsharpe[4:length(port_maxsharpe)]
}

windows()
barplot(t(wpt_roll), cex.names = 0.75,ylab="Part. (%)", 
        legend=TRUE,xlim = c(0, 130), xlab= "Portafolio óptimos")






#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestion Moderna de portafolio
# Capitulo 3: Modelo CAPM
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

# Ejemplo 3.1
## Acciones seleccionadas

fechai <- '2015-12-01'
fechaf <- '2020-12-31'
periodicidad <- "monthly"   
activo <- c("AAPL")

precios.act <- precios(activo,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.act))[-1,]
(ri = mean(retornos)*12)

indice <- c("^GSPC")
indice.hist <- precios(indice,fechai,fechaf,periodicidad)
r.indice <- diff(log(indice.hist))[-1,]
rm = mean(r.indice)
sigmam <- sd(r.indice)

datos <- coredata(cbind(retornos,r.indice))

windows()
plot(datos[,2],datos[,1],xlab="Indice S&P 500", ylab="AAPL",pch=19)
abline(lm(datos[,1]~datos[,2]), col="black") 
abline(v=0,h=0,col="gray")

# Estimacion betas

modelo<-lm(retornos~r.indice)
coef <- modelo[["coefficients"]]
(beta <- coef[2])

summary(modelo)

# ----------------------------------------------------------------

# Ejemplo 3.2

# Betas para los n activos
fechai <- '2015-12-01'
fechaf <- '2020-12-31'
periodicidad <- "monthly"   
activos <- c("AAPL","AMZN","MSFT","GOOG")

precios.act <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.act))[-1,]
mu <- colMeans(retornos)
cov <- cov(retornos)
var <- diag(cov)
sigma <- sqrt(var)

indice <- c("^GSPC")
indice.hist <- precios(indice,fechai,fechaf,periodicidad)
r.indice <- diff(log(indice.hist))[-1,]
rm = mean(r.indice)
sigmam <- sd(r.indice)

n <- length(mu)
betas <- matrix(0,ncol=n)
alphas <- matrix(0,ncol=n)
sigmaerror <- matrix(0,ncol=n)

for(i in 1:n){
    capm <- lm(retornos[,i]~r.indice)
    coef <- capm[["coefficients"]]
    betas[,i] <- coef[2]
    alphas[,i] <- coef[1]
    sigmaerror[,i] <- sd(capm$residuals)^2
}

betas
alphas

rf <- 0.0
rma <- round(rm,4)
ri <- matrix(0,ncol=n)
for(i in 1:n){
     ri[,i] <- (rf+betas[i]*(rma-rf))
}         
round(ri,4)             


#-------------------------------------------

# Ejemplo 3.3
# Modelo Sharpe - portafolio optimo

rm = mean(r.indice)
sigmam <- sd(r.indice)
sigmaeim <- t(cbind(sigmaerror,sigmam^2))
mui <- matrix(mu)

unos <- matrix(rep(1,n))
betam <- t(cbind(betas,-1))
uno <- rbind(unos,0)
alpham <-  t(cbind(alphas,rm)) 
sigmaei <- matrix(0,n+1,n+1)
diag(sigmaei) <- sigmaeim 

alb <- cbind(alpham,uno,betam)
round(alb,4)

# Matriz G
G <- t(alb)%*%solve(sigmaei)%*%alb
round(G,4)

# Vector resultado
# Ej: 
rp <- 0.02
res <- rbind(rp, 1, 0)

# Calculo pesos optimos
wpo <- solve(sigmaei)%*%alb%*%solve(G)%*%res 
wpo <- wpo[1:n]
round(wpo,4)

# Frontera eficiente
nport <- 100
j <- seq(min(mu),max(mu)+0.02, length=nport) 
wpoS <- matrix(c(0), ncol=n, nrow=nport) 
rpoS <- matrix(c(0), nrow=nport)
sigmapoS <- matrix(c(0), nrow=nport)
wj <- 0

for(i in 1:nport){
    rp <- j[i]
    res <- rbind(rp, 1, 0)
    wj <- solve(sigmaei)%*%alb%*%solve(G)%*%res 
    wj <- wj[1:n]
    wpoS[i,] <- wj
    rpoS[i,] <- wj%*%mui
    sigmapoS[i,] <- sqrt(t(wj)%*%cov%*%wj)
}

sigmaa <- c(0.0858,0.0807,0.055,0.0609)

# Plot
windows()
plot(sigmaa,mu, xlim=c(0.02,max(sigma)*1.2), ylim=c(0.01,max(mu)*1.2),cex=0.6,
     xlab="Riesgo",ylab="Retorno",font = 6)
points(sigmapo,rpo,type='l',lty=2,lwd=2)
points(sigmapoS,rpoS,type='l')
text(sigmaa,mu,labels=activos,pos = 4, cex=0.7,font = 3)
legend("bottomright",legend = c("FE Markowitz","FE modelo de mercado"),
       lty=c(2,1),cex=0.9)
grid()

# -------------------------------------------------------------

# Ejemplos 3.4 y 3.5
# Clasificacion y optimizacion Modelo de Treynor

activos <- c("AAPL","ABT","AMZN","CAT","CSX","FB","GOOG","HD",
             "JNJ","MSFT","MCD","V") 

precios.act <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.act))[-1,]

# Se toman los parámetros anualizados
mu <- colMeans(retornos)*12
cov <- cov(retornos)*12
vari <- diag(cov)
sigmai <- sqrt(vari)

rm = mean(r.indice)*12
varm <- var(r.indice)*12
sigmam <- sqrt(varm)

# Coef. Treynor
n <- length(mu)
betas <- matrix(0,ncol=n)
sigmaerror <- matrix(0,ncol=n)

for(i in 1:n){
    capm <- lm(retornos[,i]~r.indice)
    coef <- capm[["coefficients"]]
    betas[,i] <- coef[2]
    #sigmaerror[,i] <- sd(capm$residuals)^2*12 #varianzas anualizadas
}

varerror <- vari - varm%*%betas^2
treynor <- (mu-rf)/betas

matrix <- t(rbind(treynor,t(mu),t(sigmai),betas,varerror))
sort.matrix <- matrix[order(-matrix[,1]),]
ratio1 <- ((sort.matrix[,2]-rf)*sort.matrix[,4])/sort.matrix[,5]
ratio2 <- sort.matrix[,4]^2/sort.matrix[,5]
sumacu1 <- cumsum(ratio1)
sumacu2 <- cumsum(ratio2)
coef.c <- (varm%*%sumacu1)/(1+varm%*%sumacu2)

diff <- sort.matrix[,1]-coef.c
cond.diff <- diff[!is.na(diff)&diff>0]
n.optimo <- length(cond.diff)
cuttoff <- coef.c[n.optimo]

Zi <- (sort.matrix[,4]/sort.matrix[,5])*(sort.matrix[,1]-cuttoff)
Zi <- pmax(Zi,0)
(wpot <- round(Zi/sum(Zi),4))
(rpot <- t(wpot)%*%mu)
(sigmapot <- sqrt(t(wpot)%*%cov%*%wpot))

data <- cbind(sort.matrix[,1],as.vector(coef.c))
names <- rownames(data)

windows()
plot(data[,1],type="l",xlab="Número de activos", ylab="Treynor")
lines(data[,1])
points(data[,2],col="black",type="l",lty=2)
points(8,cuttoff,col="black",pch=19)
#abline(v=8,col="gray")
text(8,cuttoff,labels="C*",pos = 3, cex=1)
legend("topright",legend = c("Treynor","Tasa C"),lty=c(1,2))


matrix2 <- cbind(sort.matrix,t(coef.c))

windows()
barplot(t(wpot),beside=TRUE,col="black",xlim=c(1,15),ylab="Part. (%)")

# Comparacion Treynor-Sharpe
wsharpe <- wpt2
wtreynor <- wpot
nactivos <- names(wtreynor)

pesos <- matrix(0,2,n)
colnames(pesos) <- nactivos 
rownames(pesos) <- c("Treynor","Sharpe")
pesos[1,] <- cbind(wtreynor[nactivos[1:n]])
pesos[2,] <- cbind(wsharpe[nactivos[1:n]])

windows()
barplot(pesos,beside=TRUE,xlim=c(1,22),ylab="Part. (%)")
legend("topright",legend = c("Treynor","Sharpe"),
       fill=c("black","gray"))



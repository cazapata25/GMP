
#-----------------------------------------------------------
#-----------------------------------------------------------
# Gestión Moderna de Portafolio
# Capitulo 6: Medida Omega
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
library(ROI)

## Acciones seleccionadas
activos <- c("AAPL","AMZN","GOOG","MSFT")
fechai <- '2009-12-01'
fechaf <- '2021-12-31'
periodicidad <- "monthly"   

precios.hist <- precios(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios.hist))[-1,]

(mu <- colMeans(retornos))
cov <- cov(retornos)
var <- diag(cov)
sigma <- sqrt(var)
n <- length(mu)

##--------------------------------------------------------------
## Modelo Omega
##--------------------------------------------------------------

## FOrmulacion general: x* <- y*/z*
## Fuente: Vana, Schwendinger y Hochreiter (2016)

h <- 0

# Se define la funcion Omega

omega <- function(retornos, h){
    s <- nrow(retornos)
    Amat <- rbind(cbind(as.matrix(retornos), diag(s), 0),# u_s >= h - r_s' y
                  c(rep(0, n), rep(1, s),  0), # sum(u) = 1
                  c(rep(1, n), rep(0, s), -1), # sum(y) = z
                  c(mu,        rep(0, s), - h),
                  c(rep(0, n), rep(0, s),  1)) # mu' y  >= k * z
    var.names <- c(paste0("y_omega_aux", seq_len(n)),
                   paste0("u_omega_aux", seq_len(s)), 
                   "z_omega")
    constraint <-  L_constraint(L = Amat,  
                                dir = c(rep(">=", s), "==", "==", ">=", ">"),
                                rhs = c(rep(h, s), 1, 0, 0, 1e-05), 
                                names=var.names)
    objective  <- L_objective(L = c(mu, rep(0, s), -h))
    list(objective = objective, constraint = constraint) 
}

budget_constraint <- function(retornos, dir = "==", rhs = 1) {
    x.names <- activos
    L_constraint(L = rep(1, n), 
                 dir = dir,  rhs = rhs, names = x.names)
}

tmp <- omega(retornos,h)
pomega  <- OP(objective  =  tmp$objective,
              constraints = tmp$constraint, maximum = T)
solom <- ROI_solve(pomega, solver = "glpk")
sol_omega <- solution(solom)
wpomega <- round(sol_omega[1:n]/sol_omega["z_omega"], 4)
names(wpomega) <- activos
wpomega

windows()
barplot(t(wpomega),beside=TRUE, ylim=c(0,0.5),cex.names = 0.85,ylab="Part. (%)")

names(wpomega) <- activos
rpomega <- mu%*%wpomega
sigmapomega <- sqrt(t(wpomega)%*%cov%*%wpomega)

# Calculo medida Omega de lo activos
omegai <- NULL
for(i in 1:n){
    omegai[i] <- sum(pmax(retornos[,i]-h,0))/sum(pmax(h-retornos[,i],0))
}
names(omegai) <- activos
omegai

# --------------------------------------------------------------------------------------------------

## Formulación usando ROML
## Fuente: Vana, Schwendinger y Hochreiter (2016)

library("ROI")
library("ROML")
library("ROML.portfolio")

h <- 0
m <- model()
m$variable(portfolio, lb = 0)  
m$maximize( omega(portfolio) )
opt <- optimize(m, solver="glpk", 
                data=list(returns = coredata(retornos))) 
wpomega <- round(opt$solution[grep("portfolio", names(opt$solution))]/
                     opt$solution[grep("z", names(opt$solution))], 6)


## --------------------------------------------------------------
## Comparacion Omega y Sharpe

pesos <- cbind(wpomega,wpt)

windows()
barplot(t(pesos),beside=TRUE,ylim=c(0,0.5),cex.names = 0.85,ylab="Part. (%)")
legend("topright",legend = c("Omega","Sharpe"),
       fill=c("black","gray"),cex=0.8)


# Riesgo de los activos y del portafolio
riskomegai <- NULL
for(i in 1:n){
    riskomegai[i] <- sum(pmax(h-retornos[,i],0))
}
names(riskomegai) <- activos
riskomegai

rphist <- retornos%*%wpomega
riskomegap <- sum(pmax(h-rphist,0))

# Omega del portafolio
omegap <- sum(pmax(rphist-h,0))/sum(pmax(h-rphist,0))

## ----------------------------
# Plano Riesgo-Retorno Omega
windows()
plot(riskomegai,mu, xlim=c(min(riskomegai)*0.6,max(riskomegai)*1.1), 
     ylim=c(min(mu)*0.8,max(mu)*1.1),
     cex=0.6,xlab="Riesgo",ylab="Retorno esperado",bty="L")
text(riskomegai,mu,labels=activos,pos = 4, cex=0.7)
points(riskomegap,rpomega, col="black",pch=20)
text(riskomegap,rpomega,labels="PO",pos = 2)
abline(v=min(riskomegai),lty=2,col="gray")
arrows(min(riskomegai)*0.99, rpomega, riskomegap*1.04, rpomega)

# Efceto del umbral h:

h <- seq(0,mean(rphist),length=100)
omegap <-  matrix(0, nrow=nport)
for(i in 1:length(h)){
    m <- model()
    m$variable(portfolio, lb = 0)  
    m$maximize( omega(portfolio) )
    opt <- optimize(m, solver="glpk", 
                    data=list(returns = coredata(retornos))) 
    wpomega <- round(opt$solution[grep("portfolio", names(opt$solution))]/
                         opt$solution[grep("z", names(opt$solution))], 6)
    rphist <- retornos%*%wpomega
    omegap[i] <- sum(pmax(rphist-h[i],0))/sum(pmax(h[i]-rphist,0))
}

windows()
plot(h,omegap,xlab="Umbral h",ylab="Medida Omega",type="l",bty="L")

## ------------------------------------------------------------------




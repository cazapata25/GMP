
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Gestion Moderna de portafolio
# Funciones auxiliares
# Copyright 2022
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
# Parte 1. Modelo MV: Markowitz y Sharpe
#
#-----------------------------------------------------------------------
# 
# Optimizacion sin restricciones en los cortos

FronteraEficiente <- function(retornos){
    mu <- colMeans(retornos)
    cov <- cov(retornos)
    activos <- activos
    ones <- rep(1,n)
    x <- t(mu)%*%solve(cov)%*%mu
    y <- t(mu)%*%solve(cov)%*%ones
    z <- t(ones)%*%solve(cov)%*%ones
    d <- x*z - y*y
    g <- (solve(cov,ones)%*%x-solve(cov,mu)%*%y)%*%solve(d)
    h <- (solve(cov,mu)%*%z-solve(cov,ones)%*%y)%*%solve(d)
    rpmin = min(mu)
    rpmax <- max(mu)*1.2
    nport <- 100
    j <- seq(rpmin,rpmax, length=nport) 
    wpo <- matrix(c(0), ncol=n, nrow=nport) 
    rpo <- matrix(c(0), nrow=nport)
    sigmapo <- matrix(c(0), nrow=nport)
    wj <- 0
    cont <- 1
    for(i in 1:nport){
        wj <- g + h*j[i] 
        wpo[i,] <- t(wj)
        rpo[i,] <- t(wj)%*%mu
        sigmapo[i,] <- sqrt(t(wj)%*%cov%*%wj)
    }
    # PMVG
    wpmvg <- solve(cov,ones)%*%(1/z)
    rpmvg <- mu%*%wpmvg
    sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
    # Sharpe
    Er <- mu-rf 
    Z <- solve(cov,Er)  
    sumZ <- sum(Z) 
    wpt <- Z/sumZ 
    rpt <- t(wpt)%*%mu
    sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
    
    EF <- list()
    EF[[1]] <- wpo
    EF[[2]] <- rpo
    EF[[3]] <- sigmapo
    EF[[4]] <- cbind(wpmvg)
    EF[[5]] <- rpmvg
    EF[[6]] <- sigmapmvg
    EF[[7]] <- cbind(wpt)
    EF[[8]] <- rpt 
    EF[[9]] <- sigmapt
    return(EF)
}


 
#-----------------------------------------------------------
#-----------------------------------------------------------




#-----------------------------------------------------------

performance <- function(retornos,pesos){
    t <- nrow(retornos)
    n <- ncol(retornos)
    rport <- matrix(0,nrow=t,ncol=1)
    colnames(rport) <- c("Retornos")
    vport <- matrix(0,nrow=t,ncol=1)
    colnames(vport) <- c("Valor")
    pesos <- pesos 
    
    # Retornos historicos
    ret.hist <- retornos%*%pesos
    rport[,1] <- ret.hist
    
    # Valor del portafolio
    port.v <- matrix(0, nrow=t)
    port.v[1] <- valor
    for(i in 2:t){
        port.v[i] <- port.v[i-1]*exp(rpmv[i-1])
    }
    vport[,1] <- port.v
    
    DH <- list()
    DH[[1]] <- vport
    DH[[2]] <- rport
    return(DH)
}

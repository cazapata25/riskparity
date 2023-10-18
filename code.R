

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Paper: Paridad de riesgo y diversificacion en portafolios de inversion
# Copyright 2023
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

# Libraries
rm(list=ls())
library(zoo)
library(xts)
library(quantmod)
library(quadprog)
library(rootSolve)
library(PerformanceAnalytics)

#------------------------------------------------------------------
# Funcion para importar la informacion de Yahoo Finance
precios.hist <- function(activos,fechai,fechaf,periodicidad){
  precios <- xts()
  for(i in 1:length(activos)){
    tmp <- Ad(getSymbols(activos[i],from=fechai,to=fechaf,periodicity=periodicidad,auto.assign=FALSE))
    precios <- cbind(precios,tmp)
    precios <- na.approx(precios,na.rm=FALSE)
  }
  colnames(precios) <- activos
  tclass(precios) <- "Date"
  return(precios)
}
#------------------------------------------------------------------

## ----------------------------------------------------------------
# Lista original de activos
# activos <- c("AAPL","AMGN","BA","CAT","CSCO","CVX","DIS","HD","HON",
#             "IBM","INTC","JNJ","JPM","KO","MCD","MMM","MSFT","MRK",
#             "NKE","PG","T","TRV","UNH","V","VZ","WBA","WMT","XOM")
# Criterio de eleccion: coef. de Sharpe
# sharpei <- mu/sigma
# s <- sharpei[order(-sharpei)]
# sel <- names(s[1:15])
## -----------------------------------------------------------------

## Acciones seleccionadas
activos <- c("AAPL","AMGN","HD","HON","JNJ","KO","MCD","MRK","MSFT",
             "NKE","PG","TRV","UNH","V","WMT") 
fechai <- '2009-12-01'
fechaf <- '2022-12-31'
periodicidad <- "monthly"   
precios <- precios.hist(activos,fechai,fechaf,periodicidad)
retornos <- diff(log(precios))[-1,]

# indice <- precios.hist("^DJI",fechai,fechaf,periodicidad)
# rm <- diff(log(indice))[-1,]
# mean(rm), sd(rm)
# mean(rm)/sd(rm)

# --------------------------------------------------------------

(mu <- colMeans(retornos))
(cov <- cov(retornos))
var <- diag(cov)
(sigma <- sqrt(var))
n <- length(mu)
(sharpei <- round(mu/sigma,4))

# ------------------------------------------------------------------
# Comparacion de resultados

# Solucion de los portafolios optimos: MV y PR

# PR Naive
(wrpnaive <- 1/sigma/(sum(1/sigma)))

# PR vanilla
b <- rep(1/n,n)
f_root <- function(x,Sigma){
  cov <- Sigma
  n <- ncol(cov)
  return(cov%*%x - b/x)
}

x_root <- multiroot(f_root, start=b, parms=cov)$root
wRPv <- x_root/sum(x_root)
names(wRPv) <- activos
wRPv

# MV de Markowitz
unos <- rep(1, n)
solMVG <- solve.QP(Dmat=cov*2, dvec=rep(0,n), Amat=cbind(unos,diag(1,n)),
                   bvec=c(1,rep(0,n)), meq=1)
wpmvg <- round((solMVG$solution/sum(solMVG$solution)),4)
names(wpmvg) <- activos
wpmvg

# Figura 1
pesos <- rbind(wRPv,wrpnaive,wpmvg)

windows()
barplot(pesos,beside=TRUE,ylim=c(0,0.3),ylab="Part. (%)") #
legend("topleft",c("PRn","PRv","MV"), fill=c("black","darkgray","gray"),bty="n")

# Graficos de desempeno

valor <- 1 
t <- length(rm)
#MV
rpmv <- c(retornos%*%wpmvg)
#RP naive
rprp <- c(retornos%*%wrpnaive)
#RP vanilla
rprpv <- c(retornos%*%wRPv)
# Retornos acumulados
racummv<- rep(0,t)  
racumrp <- rep(0,t)
racumrpv <- rep(0,t)
racumindice <- rep(0,t)
racummv[1] <- valor
racumrp[1] <- valor
racumrpv[1] <- valor
racumindice[1] <- valor

for(i in 2:t){
  racummv[i] <- racummv[i-1] * exp(rpmv[i-1])
  racumrp[i] <- racumrp[i-1] * exp(rprp[i-1])
  racumrpv[i] <- racumrpv[i-1] * exp(rprpv[i-1])
  racumindice[i] <-  racumindice[i-1] * exp(rm[i-1])
}

# Figura 2

windows()
plot(racumrp, col="black",type="l", ylim=c(0.8,max(racummv,racumrp)),
     ylab="Retorno acumulado",xlab="Tiempo (meses)")
lines(racumrpv, col="black",lty=2)
lines(racummv, col="darkgray")
lines(racumindice, col="darkgray",lty=2)
legend("topleft",c("PRn","PRv","MV","DJI"),
       lty =c(1,2,1,2), col=c("black","black","darkgray","darkgray"),bty="n")

# --------------------------------------------------------

# Consolidado de los resultados
rp <- c(mean(rprp),mean(rprpv),mean(rpmv),mean(rm))
sd <- c(sd(rprp),sd(rprpv),sd(rpmv),sd(rm))
sharp <- c((mean(rprp)/sd(rprp)),(mean(rprpv)/sd(rprpv)),(mean(rpmv)/sd(rpmv)),(mean(rm)/sd(rm)))
nactivos <- c(length(wrpnaive[wrpnaive>0]),length(wRPv[wRPv>0]),length(wpmvg[wpmvg>0]),29)
resumen <- rbind(round(rp,4),round(sd,4),round(sharp,4),nactivos)
rownames(resumen) <- c("Retorno","Volatilidad","Coef. Sharpe","No. Activos")
colnames(resumen) <- c("PRn","PRv","MV","DJI")
resumen

# -------------------------------------------------------------
# -------------------------------------------------------------

# Rolling-window analisis 

t <- 60
tf <- nrow(retornos)-12
i <- c(t:tf)
wprp_roll <- matrix(0,nrow=length(i),ncol=n)
wpmv_roll <- matrix(0,nrow=length(i),ncol=n)
colnames(wprp_roll) <- activos
colnames(wpmv_roll) <- activos
rpmv <- rep(0,length(i))
rprp <- rep(0,length(i))
rb <- rep(0,length(i))
sdmv <- rep(0,length(i))
sdrp <- rep(0,length(i))
sdrb <- rep(0,length(i))
sharpemv <- rep(0,length(i))
sharperp <- rep(0,length(i))
sharperb <- rep(0,length(i))

for(j in 1:length(i)){
  ret_roll <- retornos[1:i[j],] 
  mu_roll <- colMeans(ret_roll)
  cov_roll <- cov(ret_roll)
  sigma_roll <- sqrt(diag(cov_roll))
  rpm <- rm[1:i[j],]
  wrp <- 1/sigma_roll/(sum(1/sigma_roll))
  solMVG <- solve.QP(Dmat=cov_roll*2, dvec=rep(0,n), Amat=cbind(unos,diag(1,n)),
                     bvec=c(1,rep(0,n)), meq=1)
  wpmv <- round((solMVG$solution/sum(solMVG$solution)),6)
  wprp_roll[j,] <- rbind(wrp)
  wpmv_roll[j,] <- rbind(wpmv)
  ret_out <- ret_roll 
  rpmv[j] <- mean(ret_roll%*%wpmv) #MV
  rprp[j] <- mean(ret_roll%*%wrp) #RP
  rb[j] <- mean(rpm)    #DJI
  sdmv[j] <- sd(ret_roll%*%wpmv)
  sdrp[j] <- sd(ret_roll%*%wrp)
  sdrb[j] <- sd(rpm) 
  sharpemv[j] <- mean(ret_roll%*%wpmv)/sd(ret_roll%*%wpmv)
  sharperp[j] <- mean(ret_roll%*%wrp)/sd(ret_roll%*%wrp)
  sharperb[j] <- mean(rpm)/sd(rpm)
}

windows()
barplot(t(wprp_roll), cex.names = 0.75,ylab="Part. (%)", 
        legend=TRUE,xlim = c(0, 130), xlab= "Portafolio optimos")

windows()
barplot(t(wpmv_roll), cex.names = 0.75,ylab="Part. (%)", 
        legend=TRUE,xlim = c(0, 130), xlab= "Portafolio optimos")

windows()
plot(rprp,ylim=c(0.004,max(rprp,rpmv)),type="l",ylab="Retorno esperado",
     xlab="Rebalanceo de los portafolios")
lines(rpmv,col="black",lty=2)
lines(rb,col="gray")
legend("bottomright",c("PRn","MV","DJI"),
       lty =c(1,2,1), col=c("black","black","darkgray"),bty="n")

windows()
plot(sharpemv,ylim=c(0.1,max(sharperp,sharpemv)),type="l",
     ylab="Coeficiente de Sharpe",xlab="Rebalanceo de los portafolios")
lines(sharperp,col="black",lty=2)
lines(sharperb,col="gray")
legend("bottomright",c("PRn","MV","DJI"),
       lty =c(1,2,1), col=c("black","black","darkgray"),bty="n")

# ---------------------------------------------------------------------
# HH Index

hhiMV <- matrix(0,nrow=nrow(wpmv_roll))
for(i in 1:nrow(wpmv_roll)){hhiMV[i] <- sum((wpmv_roll[i,]*100)^2)}

hhiPR <- matrix(0,nrow=nrow(wprp_roll))
for(i in 1:nrow(wprp_roll)){hhiPR[i] <- sum((wprp_roll[i,]*100)^2)}

windows()
plot(hhiMV,xlab="Portafolios", ylab="HHI", type="l", lty=2, 
     ylim=c(0,2500))
lines(hhiPR)
legend("topright",c("MV","PRn"),lty =c(2,1), 
       lwd =c(2,1), col=c("black","black"),bty="n")


# --------------------------------------------------------------
# End ...
# --------------------------------------------------------------

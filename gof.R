
rm(list=ls())

library(goftest)
library(dplyr)
library(tidyr)

#=======================================================
findCV <- function(tvec,tau,sigma,cv){
        n=length(tvec)
        xvec=diff(c(0,tvec))
        if(sigma=="s"){
                sdest <- sqrt(var(xvec))
                mu <- mean(xvec)
                CV <- sdest/mu
        }
        if(sigma=="c"){
                xvec2 <- c(xvec,tau-tvec[n])
                mu <- tau/n
                sdest <- sqrt(sum(xvec2^2)/n-mu^2)
                CV <- sdest/mu
        }
        if(sigma=="l"){
                xdiff <- xvec[2:n]-xvec[1:(n-1)]
                sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
                mu <- mean(xvec)
                CV <- sdest/mu
        }
        if(sigma=="fixCV"){
                CV <- cv
        }
        return(CV)
}  

## Lewis Robinson test

LRtestobs <- function(tvec,tau,sigma,cv){
        CV <- findCV(tvec,tau,sigma,cv)
        n=length(tvec)
        LR <- (sqrt(12)/(CV*tau*sqrt(n)))*(sum(tvec)-n*tau/2)
        return(LR)
}


LRtest <- function(tvec,tau,sigma="s",cv=1){
        LR <- LRtestobs(tvec,tau,sigma,cv)
        pLR <- 2*pnorm(-abs(LR))
        cat("Lewis-Robinson test for trend in time censored processes.\n")
        cat(paste("Test statistic: LR =",round(LR,digits=3),"\n"))
        cat(paste("p-value =",round(pLR,digits=5)),"\n")
        CV <- findCV(tvec,tau,sigma,cv)
        cat(paste("CV =", round(CV,digits=3) ,"\n"))
}


## Cramer von-Mises test

CvMtestobs <- function(tvec,tau,sigma,cv){
        CV <- findCV(tvec,tau,sigma,cv)
        n=length(tvec)
        xvec=diff(c(0,tvec))
        konst <- 1/(CV^2*n)
        indseq <- 0:(n-1)
        sumledd <- sum(indseq^2*xvec/tau)-n*sum(indseq*(tvec^2-c(0,tvec[1:(n-1)])^2)/tau^2)
        sisteledd <- n^2/3+n^2*(tvec[n]^2/tau^2-tvec[n]/tau)
        CV <- konst*(sumledd+sisteledd)
        return(CV)
}  


CvMtest <- function(tvec,tau,sigma="s",cv=1){
        CvM <- CvMtestobs(tvec,tau,sigma,cv)
        pCvM <- 1-pCvM(q=CvM)
        cat("Cramer von-Mises test for trend in time censored processes.\n")
        cat(paste("Test statistic: CvM =",round(CvM,digits=3),"\n"))
        cat(paste("p-value =",round(pCvM,digits=5)),"\n")
        CV <- findCV(tvec,tau,sigma,cv)
        cat(paste("CV =", round(CV,digits=3) ,"\n"))
}

#=======================================================



d <- read.csv("cd_cars.csv")
d <- as_tibble(d)
dd <- d %>% arrange(FM, Mil)
LHD_tau <- 3000


#=====================================================
# Nelson-Aalen plot for trend:
#=====================================================

LHD1_ftimes <- dd$Mil[1:76]
LHD1_N <- length(LHD1_ftimes) 
plot(LHD1_ftimes,1:LHD1_N,xlab="Mileage",ylab="Cumulative number of failures")
lines(c(0,LHD1_ftimes[LHD1_N]),c(0,LHD1_N),lty=2)


LHD2_ftimes <- dd$Mil[77:163]
LHD2_N <- length(LHD2_ftimes) 
plot(LHD2_ftimes,1:LHD2_N,xlab="Mileage",ylab="Cumulative number of failures")
lines(c(0,LHD2_ftimes[LHD2_N]),c(0,LHD2_N),lty=2)


LHD3_ftimes <- dd$Mil[164:274]
LHD3_N <- length(LHD3_ftimes) 
plot(LHD3_ftimes,1:LHD3_N,xlab="Mileage",ylab="Cumulative number of failures")
lines(c(0,LHD3_ftimes[LHD3_N]),c(0,LHD3_N),lty=2)


#=========================================
# Run the trend tests
# Warranty claim dataset (cd_cars.csv)
#=========================================

LRtest(LHD1_ftimes,LHD_tau,sigma="s")
LRtest(LHD1_ftimes,LHD_tau,sigma="l")

CvMtest(LHD1_ftimes,LHD_tau,sigma="s")
CvMtest(LHD1_ftimes,LHD_tau,sigma="l")


LRtest(LHD2_ftimes,LHD_tau,sigma="s")
LRtest(LHD2_ftimes,LHD_tau,sigma="l")

CvMtest(LHD2_ftimes,LHD_tau,sigma="s")
CvMtest(LHD2_ftimes,LHD_tau,sigma="l")


LRtest(LHD3_ftimes,LHD_tau,sigma="s")
LRtest(LHD3_ftimes,LHD_tau,sigma="l")

CvMtest(LHD3_ftimes,LHD_tau,sigma="s")
CvMtest(LHD3_ftimes,LHD_tau,sigma="l")

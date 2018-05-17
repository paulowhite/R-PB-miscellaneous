### Generate-toy-data-for-permutation-test.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: May 16 2018 (09:07) 
## Version: 
## Last-Updated: May 17 2018 (09:26) 
##           By: Paul Blanche
##     Update #: 56
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

rm(list=ls())

library(tidyr)

set.sed(4658)

## {{{ parameter for data generation
n1 <- n0 <- 10     # size of the 2 groups
time <- c(-1:2,6:8) # time of the measurements
mysd <- 1.5 # sd of the noise
mysdRandomInt <- 1.5 # ds of random intercept
propmissing <- 0.25  # proportion of missing marker measure values
f1 <- function(x){1-0.5*x}   # function for mean trjectory for treated
f0 <- function(x){0+0.35*x}  # function for mean trjectory for treated
## }}}


## {{{ create data
groupwide <- c(rep(1,n1),rep(0,n0)) # 
d <- data.frame(id=rep(1:(n1+n0),each=length(time)),
                group=rep(groupwide,each=length(time)),
                time=rep(time,times=(n0+n1)))
d$Y <- do.call("c",lapply(groupwide,function(x){x*f1(time) + (1-x)*f0(time) + rnorm(length(time),sd=mysd)  }))
d$Y <- d$Y+ rep(rnorm(n=(n1+n0),mean=0,sd=mysdRandomInt),each=length(time)) # add random effect
# create unbalanced data, by deleteing some lines (about propmissing%)
d <- d[-sample(x=1:nrow(d),size=floor(propmissing*nrow(d)),replace=FALSE),]
d0 <- d[d$group==0,]
d1 <- d[d$group==1,]
dnew <- d
## }}}

## {{{ plot the generated data
interaction.plot(d0$time,
                 d0$id,
                 d0$Y,
                 xlab="time",
                 ylab="Y",
                 legend=F,col="blue",
                 ylim=range(d$Y))
interaction.plot(d1$time,
                 d1$id,
                 d1$Y,
                 xlab="time",
                 ylab="Y",
                 legend=F,
                 add=TRUE,
                 col="red")
Mean0 <- tapply(d0$Y,d0$time,mean)
lines(x=1:length(Mean0),Mean0, type="b",col="blue",lwd=2)
Mean1 <- tapply(d1$Y,d1$time,mean)
lines(x=1:length(Mean1),Mean1, type="b",col="red",lwd=2)
legend("topleft",legend=paste0("group=",0:1),col=c("blue","red"),lwd=1,bty="n")
## }}}

## {{{ rename with stuid names (to check the robustness of the function to variable names)
names(dnew) <- c("MYID","MYGROUP","THETIME","THEY")
## }}}


## {{{ Save toy data
write.table(dnew,
            "~/research/Rigshospitalet/Stefan/Pig-longitudinal-trajectories/test-trajectories/toy-data-1.csv",
            row.names=FALSE,
            sep=";")
## ddd <- read.table("~/research/Rigshospitalet/Stefan/Pig-longitudinal-trajectories/test-trajectories/toy-data-1.csv",
                  ## sep=";",
                  ## header=TRUE)
## head(ddd)
## }}}



##----------------------------------------------------------------------
### Generate-toy-data-for-permutation-test.R ends here

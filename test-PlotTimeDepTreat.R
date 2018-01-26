### test-PlotTimeDepTreat.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Jan 26 2018 (14:24) 
## Version: 
## Last-Updated: Jan 26 2018 (14:28) 
##           By: Paul Blanche
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

rm(list=ls())

setwd("~/research/Gentofte/Frederik/Time-dep-treatment/")

source("PlotTimeDepTreat.R")

## {{{ load data
ddd <- read.table(file="d1.csv",
                sep=";", header = TRUE)
head(ddd)
runif(1)
names(ddd)
summary(ddd)
## }}}

## {{{ old test
## PlotTimeDepTreat(data=ddd,                   # data to plot
## N=20,                   # subsample size
## seed=1234,              # random seed to reproducibility (random subsample))
## name.id="pnr",           # name of id variable´
## name.time="time",         # name of time to event variable
## name.status="status",       # name of status variable
## name.dateEntry="dateEntry",    # name of subject enrollment date variable
## name.treat="Treat",        # name treatment variable
## name.start="starts",        # name variable of the date when treatment starts
## name.stop="stops",         # name variable of the date when treatment stops
## ##                 treat.levels,      # names of all the levels of the treatment variable (to assign the colors)
## status.levels=c("censored","death","event"),     # names of all the levels of the status variable (to assign the symbols)
## pch.status=c(1,19,17),        # pch to display event types
## cex.status=2,           # cex to display event types
## lwd=3,                  # thickness of the horizontal lines
## ##            col,               # colors to display treatment levels
## add.warning=TRUE        # default is to warn that this displays individual data and that therefore we cannot take this plot out of  Denmark Statistics
## )

## }}}

## {{{ plot
x <- PlotTimeDepTreat(data=ddd,                   # data to plot
                      N=50,                   # subsample size
                      seed=12345,              # random seed to reproducibility (random subsample))
                      name.id="id",           # name of id variable´
                      name.time="T",         # name of time to event variable
                      name.status="delta",       # name of status variable
                      name.dateEntry="fwstart",    # name of subject enrollment date variable
                      name.treat="exposure",        # name treatment variable
                      name.start="begins",        # name variable of the date when treatment starts
                      name.stop="finishes",         # name variable of the date when treatment stops
                      treat.levels=paste("Treat",1:5,sep="."),      # names of all the levels of the treatment variable (to assign the colors)
                      status.levels=c("censored","death","event"),     # names of all the levels of the status variable (to assign the symbols)
                      pch.status=c(1,19,17),        # pch to display event types
                      cex.status=2,           # cex to display event types
                      lwd=3,                  # thickness of the horizontal lines
                      col=c("black","red","blue","orange","ForestGreen"),               # colors to display treatment levels
                      add.warning=!TRUE        # default is to warn that this displays individual data and that therefore we cannot take this plot out of  Denmark Statistics
                      )
head(x)
tail(x)
## }}}

##----------------------------------------------------------------------
### test-PlotTimeDepTreat.R ends here

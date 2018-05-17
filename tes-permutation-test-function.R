### tes-permutation-test-function.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: May 16 2018 (09:10) 
## Version: 
## Last-Updated: May 17 2018 (09:27) 
##           By: Paul Blanche
##     Update #: 23
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

rm(list=ls())


setwd("~/research/Rigshospitalet/Stefan/Pig-longitudinal-trajectories/test-trajectories/")

## {{{ load permutation test function
source("permutation-test-function.R")
## }}}


## {{{ load toy data
d <- read.table("toy-data-1.csv",
                sep=";",
                header=TRUE)
head(d)
table(d$time)
## }}}

## {{{ run function
res <- PermutTestCompare2Curves(data=d,
                                NMC=5000,
                                VarY="THEY",
                                VarGroup="MYGROUP",
                                VarID="MYID",
                                varTime="THETIME")#,
                                ## times=c(-1,1,2,8))

res$CompTime
res$NMC
res$pval
res$error
res$CIp
## }}}





##----------------------------------------------------------------------
### tes-permutation-test-function.R ends here

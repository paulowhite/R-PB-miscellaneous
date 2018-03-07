### ComparePCA.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar  7 2018 (13:34) 
## Version: 
## Last-Updated: Mar  7 2018 (14:37) 
##           By: Paul Blanche
##     Update #: 82
##----------------------------------------------------------------------
## 
### Commentary: 
##
## Aim is to compare Percent Annual Change between two groups
##
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

# Input:
# x, y        : objects fitted from LogLinDSRTrend
# conf.level  : confidence level of confidence intervals
# NMC=50000   : number of Monte Carlo simulation to compute Melted intervals
## seed=1234  : seed for reproducibility of melted confidence intervals
#
# Output:
# pval.beta : p-value for H0: PCA model x = PCA model y, based on normal approximation on beta coef of log linear model
# pval.PCA  : p-val for same H0, using delta-method, based on normal approximation on PCA
# pval.melted : p-val for same H0, based on melted confidence interval (experimental...)
# difference  : difference in PCA (between the two models)
# CI          : confidence interval for the difference in PCA, based on delta-method and on normal approximation on PCA
# CI.melted   : melted confidence interval (based on normal approx on beta of the log linear models) (experimental...)
# PCA         : PCA of the two model
# CI.Wald.PCA.x : confint PCx, normal approx on PCA scale 
# CI.Wald.PCA.y : confint PCy, normal approx on PCA scale
# CI.PCAx       : confint PCx, normal approx on beta of log linear model x
# CI.PCAy       : confint PCy, normal approx on beta of log linear model y


ComparePCA <- function(x,y,conf.level=0.95, NMC=50000, seed=1234){

    set.seed(seed) # random seed for MonteCarlo for melted  confidence interval

    ## {{{ Main stuff 
    # type I error
    alpha <- (1-conf.level)
    # beta of log.lm model
    betax <- coef(x$log.lm.fit)[2]
    betay <- coef(y$log.lm.fit)[2]
    # convert into PCA
    PCAx <- (exp(betax)-1)*100
    PCAy <- (exp(betay)-1)*100
    # s.e. of beta
    sebetax <- summary(x$log.lm.fit)$coef[2,"Std. Error"]
    sebetay <- summary(y$log.lm.fit)$coef[2,"Std. Error"]
    # s.e. of PCA via delta method
    sePCAx <- 100*exp(betax)*sebetax
    sePCAy <- 100*exp(betay)*sebetay
    # difference in PCA and s.e. (simple two sample cases)
    diff <-  PCAx - PCAy   
    sediff <- sqrt(sePCAx^2 + sePCAy^2)
    # Wald confint
    lower <- diff - qnorm(1-alpha/2)*sediff
    upper <- diff + qnorm(1-alpha/2)*sediff
    CIPCAx <- c(PCAx - qnorm(1-alpha/2)*sePCAx,
                PCAx + qnorm(1-alpha/2)*sePCAx)
    CIPCAy <- c(PCAy - qnorm(1-alpha/2)*sePCAy,
                PCAy + qnorm(1-alpha/2)*sePCAy)
    # pval on PCA or beta scales
    pvalPCA <- min(1,2*(1-pnorm(abs(diff)/sediff)))
    pvalbeta <- min(1,2*(1-pnorm(abs(betax-betay)/sqrt(sebetax^2 + sebetay^2))))
    # prepare output
    myPCA <- c(PCAx,PCAy)
    names(myPCA) <- c("x","y")
    myCI <- c(lower,upper)
    names(myPCA) <- paste(c(alpha/2,1-alpha/2)*100,"%",sep="")
    ## }}}
    
    ## {{{ Wald on beta back transformed 
    CIx <- function(alpha){
        exp(c(betax-qnorm(1-alpha/2)*sebetax,
              betax+qnorm(1-alpha/2)*sebetax))-1
        ## (exp(confint.default(x$log.lm.fit,level= 1-alpha)[2,])-1)
    }
    CIy <- function(alpha){
        exp(c(betay-qnorm(1-alpha/2)*sebetay,
              betay+qnorm(1-alpha/2)*sebetay))-1
        ## (exp(confint.default(y$log.lm.fit,level= 1-alpha)[2,])-1)
    }
    ## browser()
    ## }}}
    
    ## {{{ melted confidence interval for the difference
    ## browser()

    # melted interval, see :
    ## Fay, Michael P., Michael A. Proschan, and Erica Brittain. "Combining
    ## one‐sample confidence procedures for inference in the two‐sample case."
    ## Biometrics 71.1 (2015): 146-156.

    ForMelted <- function(alpha){
        CI <- rbind(CIy(alpha),CIx(alpha))
        colnames(CI) <- c("lower","upper")
        rownames(CI) <- c("y","x")
        CI
    }

    allsimu <- lapply(1:NMC,function(x){ForMelted(runif(1))})
    ## allsimu[[27498]]
    ## allsimu
    alllci0 <- unlist(lapply(allsimu,"[",1))
    alluci0 <- unlist(lapply(allsimu,"[",3))
    alllci1 <- unlist(lapply(allsimu,"[",2))
    alluci1 <- unlist(lapply(allsimu,"[",4))

    gGreater <- alllci1-alluci0
    summary(gGreater)
    gLess <- alluci1-alllci0   
    summary(gLess)
    ## browser()
    pgr <- length(gGreater[gGreater <= 0])/NMC
    pless <- length(gLess[gLess >= 0])/NMC
    pvalueMelted <- min(1, 2 * pgr, 2 * pless)

    meltedCI <- c(quantile(gGreater, probs = alpha),
                  quantile(gLess, probs = 1 - alpha))
    names(meltedCI) <- c("lower","upper")

    ## }}}
    
    ## {{{ output
    out <- list(pval.beta=pvalbeta,
                pval.PCA=pvalPCA,
                pval.melted=pvalueMelted,
                difference=diff,
                CI=myCI,
                CI.melted=meltedCI*100,
                PCA=myPCA,
                CI.Wald.PCA.x=CIPCAx,
                CI.Wald.PCA.y=CIPCAy,
                CI.PCAx=CIx(alpha)*100,
                CI.PCAy=CIy(alpha)*100
                )
    ## }}}    
    out
}

##----------------------------------------------------------------------
### ComparePCA.R ends here

### Compute-log-lin-trend-for-DSR.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar  6 2018 (10:04) 
## Version: 
## Last-Updated: Mar  7 2018 (14:41) 
##           By: Paul Blanche
##     Update #: 86
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## Reference: Section 4 of Fay et al., Biometrics 2006
## (Estimating average annual percent change for disease rates
## without assuming constant change).


# Input:
#
# Count      : Matrix of counts, columns=years, rows= groups (to standardize). Should have years as colnames!
# Pop        : Matrix of population, columns=years, rows= groups (to standardize). Should have years as colnames!
# ssdpop     : vector of standard population (per groups to standardize). Should have years as colnames!
# conf.level : confidence level of confidence intervals
#
#
# Output:
#
# DSR         : (Yearly) Directly standardized Rates
# CI.DSR      : Their corresponding gamma intervals (Fay et al. 1997)
# log.lm.fit  : fit of the log linear model
# PCA         : Percent Change annualized (See Section 4 of Fay et al., Biometrics 2006)
# CI.PCA      : confidence intervals for PCA (normal approx on beta of log linear model)


LogLinDSRTrend <- function(Count,Pop,sdpop,conf.level=0.95){
    ## browser()

    ## {{{ compute  key quantities
    si <- sdpop/sum(sdpop) # standard for age-group i standardized so that \sum_i si = 1
    yij <- Count           # counts for age group i (row) for year j (column)
    nij <- Pop #  mid-year population for age group i (row) for year j (column)
    lambdaij <- Count/Pop # (crude) rate for  group i (row) for year j (column)
    DSRj <- t(t(lambdaij)%*%si)  # d_j=\sum_i s_i lambdaij, just to be sure that it equals DSR, the Directly standardized Rate for year j
    varDSRj <-  t(t((yij/(nij^2)))%*%(si^2))  # estimator of the variance, defined as vj = \sum_i (si/ nij)^2 yij = \sum_i si^2 (yij / (nij)^2)    of the Directly standardized Rate for year j, that is d_j=\sum_i s_i lambdaij
    varepsilonj <-  as.vector(varDSRj/(DSRj^2)) #sigmaj^2=vj/dj^2, estimator of var(epsilonj) from delta method
    zj <- as.vector(log(DSRj)) # zj=log(dj)) 
    tj <- as.numeric(0:(length(DSRj)-1)) # time since first year
    ## }}}
    ## {{{ linear model and PCA
    fit.wls <- lm(zj~tj, weights=1/varepsilonj)
    PCA <- (exp(coef(fit.wls)[2])-1)*100 # linear mododel estimator (LME) of the PCA (percent change annualized)
    CIPCA <- (exp(confint.default(fit.wls,level= conf.level)[2,])-1)*100 # confint for PCA. Be careful, confint() whould call onfint-lm, for which the degree of freedom seem wrong. confint.default() use normal distribution instead of student, which seems more appropriate for our aggregated data. 
    predDSR <-  exp(predict(fit.wls))
    ## }}}
    ## {{{ Confidence interval for DSR

    # !!!!!!!!!! NEED TO ACCOUNT FOR MANY DSR !!!!!!
    alpha <- (1-conf.level)
    wMj <- apply(matrix(si,nrow=nrow(nij),ncol=ncol(nij),byrow = FALSE)/nij,2,max)
    lci <- qgamma(alpha/2,
                  shape = (DSRj^2)/varDSRj,
                  scale = varDSRj/DSRj)
    uci <- qgamma(1 - alpha/2,
                  shape = ((DSRj + wMj)^2)/(varDSRj + wMj^2),
                  scale = (varDSRj + wMj^2)/(DSRj + wMj))
    ## rbind(lci,DSRj,uci)
    CIDSRj <- rbind(lci,uci)
    colnames(CIDSRj) <- colnames(DSRj)
    rownames(CIDSRj) <- paste(c(alpha/2,1-alpha/2)*100,"%",sep="")
    ## }}}
    ## {{{ output
    DSRjout <- as.vector(DSRj)
    names(DSRjout) <- colnames(DSRj)
    out <- list(DSR = DSRjout,
                CI.DSR = CIDSRj,
                log.lm.fit = fit.wls,
                PCA = PCA,
                CI.PCA = CIPCA
                )
    ## }}}
    out
}

##----------------------------------------------------------------------
### Compute-log-lin-trend-for-DSR.R ends here

### plotLogLinTredforDSR.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar  7 2018 (11:13) 
## Version: 
## Last-Updated: Mar  7 2018 (14:56) 
##           By: Paul Blanche
##     Update #: 43
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


plotLogLinTredforDSR <- function(x, multiplier,
                                 add=FALSE, col="black",
                                 lty=1,
                                 lwd=1,
                                 ylim,
                                 xlim,
                                 printPCA=TRUE,
                                 digitsPCA=1,
                                 yPCA){


    if(missing(multiplier)){multiplier <- 1}
    if(missing(ylim)){ ylim <- range(x$CI.DSR*multiplier)}
    if(missing(xlim)){ xlim <- range(as.numeric(names(respsyk$DSR)))}
    if(missing(yPCA)){yPCA <- 0.95*max(ylim)}
    
    if(multiplier==1){
        myylab <- "Standadized Incidence Rate (/person-year)"
    }else{
        myylab <- paste0("Standadized Incidence Rate (/",
                         format(multiplier,scientific=TRUE)," person-years)")
                         ## format(multiplier,scientific=FALSE)," person-years)")
    }
    
    if(!add){
        plot(NULL,
             xlim=xlim,
             ylim=ylim,
             ylab=myylab,
             xlab="Year",
             axes=FALSE)
    
        axis(1,at=min(xlim):max(xlim))
        axis(2,las=2)
    }
    for(i in min(xlim):max(xlim)){
        col1 <- which(i==min(xlim):max(xlim))
        points(x=i,
               y=x$DSR[col1]*multiplier,
               pch=20)   
        arrows(angle=90,
               length=0.05,
               code=3,
               x0=i,
               x1=i,
               y0=x$CI.DSR[1,col1]*multiplier,
               y1=x$CI.DSR[2,col1]*multiplier
               )
    }
    ## browser()
    predSDR <- exp(predict(x$log.lm.fit))*multiplier
    lines(min(xlim)+0:(length(x$DSR)-1),
          predSDR,
          type="b",lwd=lwd, col=col, lty=lty)
    if(printPCA){       
        text(x=min(xlim),y=yPCA,
             paste("PCA = ",
                   round(x$PCA,digits=digitsPCA),"%,  95% CI=[",
                   round(x$CI.PCA[1],digits=digitsPCA),
                   " ; ",
                   round(x$CI.PCA[2],digits=digitsPCA),
                   "].",
                   sep=""),
             pos=4, col=col)
    }
}

##----------------------------------------------------------------------
### plotLogLinTredforDSR.R ends here

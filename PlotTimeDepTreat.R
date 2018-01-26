### PlotTimeDepTreat.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Jan 24 2018 (09:42) 
## Version: 
## Last-Updated: Jan 26 2018 (14:29) 
##           By: Paul Blanche
##     Update #: 166
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

PlotTimeDepTreat <- function(data,                   # data to plot
                             N=30,                   # subsample size
                             seed=1234,              # random seed to reproducibility (random subsample))
                             name.id,           # name of id variable´
                             name.time,         # name of time to event variable
                             name.status,       # name of status variable
                             name.dateEntry,    # name of subject enrollment date variable
                             name.treat,        # name treatment variable
                             name.start,        # name variable of the date when treatment starts
                             name.stop,         # name variable of the date when treatment stops
                             treat.levels,      # names of all the levels of the treatment variable (to assign the colors)
                             status.levels,     # names of all the levels of the status variable (to assign the symbols)
                             pch.status,        # pch to display event types
                             cex.status=2,           # cex to display event types
                             lwd=3,                  # thickness of the horizontal lines
                             col,               # colors to display treatment levels
                             add.warning=TRUE        # default is to warn that this displays individual data and that therefore we cannot take this plot out of  Denmark Statistics
                             ){

    ## browser()
    require(plyr)
    old <- .Random.seed
    set.seed(seed)

    ## browser()
    
    ## {{{ data preparation
    allthepnr <- unique(data[,name.id])
    thepnr <- allthepnr[sample(x=1:length(allthepnr),size=N)] # randomly select N subjects
    d <- data[data[,name.id] %in% thepnr,]


    AllEntries <- rep(NA,N)
    ## browser()
    d$date.Event <- as.Date(d[,name.dateEntry]) + d[,name.time]*365.25 # define date of event
    ## d[,c(name.dateEntry,"date.Event",name.time)]   
    Startfw <- min(as.Date(d[,name.dateEntry])) # start of follow-up
    Stopfw <- max(as.Date(d[,"date.Event"]))   # end of follow-up
    StudyLength <- as.numeric(difftime(Stopfw,Startfw,unit="days"))/365.25
    ## browser()
    for(i in 1:N){
        AllEntries[i] <- unique(as.character(d[d[,name.id]==thepnr[i], name.dateEntry]))
    }
    ## browser()
    timeSinceStart <- as.numeric(difftime(as.Date(AllEntries),Startfw,unit="days"))/365.25
    tosort <- order(timeSinceStart)
    thepnr <- thepnr[tosort]
    timeSinceStart <- timeSinceStart[tosort]

    if(missing(col)){
        mycols <- 1:length(unique(d[,name.treat]))
    }else{
        mycols <- col    
    }


    
    if(missing(status.levels)){status.levels <- unique(d[,name.status])}
    if(missing(pch.status)){pch.status <- 1:length(status.levels)}
    if(missing(treat.levels)){treat.levels <- unique(d[,name.treat])}

    ## browser()
    d$pchstatus <- plyr::mapvalues(d[,name.status], from = status.levels, to = pch.status, warn_missing = FALSE)
    d$mycol <- plyr::mapvalues(d[,name.treat], from = treat.levels, to = col, warn_missing = FALSE)
    ## }}}

    ## {{{ Initialize plot 
    ## browser()
    ## max(d[,name.time])
    plot(timeSinceStart,
         1:N,
         xlab="Calendar time",
         xlim=c(0,StudyLength),
         ylab="",
         axes=FALSE,
         type='n'
         )
    axis(1,at=c(0,StudyLength),labels=c(Startfw,Stopfw))   
    ## }}}

    ## {{{ loop over the subjects
    for(i in 1:N){
        cat(paste("Draw for subject",i,"out of",N,"\n"))
        ## if(i==10){ browser()}       
        di <- d[which(d[,name.id]==thepnr[i]),]
        nn <- nrow(di)
        timeEvent <- timeSinceStart[i] + di[1,name.time]
        ## browser()
        for(j in 1:nrow(di)){
            t0 <- as.numeric(difftime(as.Date(di[j,name.start]),Startfw,unit="days"))/365.25 
            t1 <- as.numeric(difftime(as.Date(di[j,name.stop]),Startfw,unit="days"))/365.25
            if(j==1){text(x=t0,y=i,labels=paste("(",nrow(di),")",sep=""),pos=2)}
            cat(paste("Period",j,"col=",di$mycol[j],"\n"))
            ## if(i*j==1){browser()}
            arrows(x0=t0,
                   x1=t1,
                   y0=i,
                   y1=i,
                   col=as.character(di$mycol[j]),
                   angle=90,
                   length=0.05,
                   code=3,
                   lwd=lwd)        
        }
        points(x=timeEvent,
               y=i,
               pch=as.numeric(as.character(di$pchstatus[1])),
               cex=cex.status)
    }


    ## }}}

    ## {{{ add legend
    legend("topleft",
           col=mycols,
           lty=1,ncol=1,
           lwd=lwd,
           legend=treat.levels,
           title="Treatment",
           bty='n')


    legend("left",
           cex=1,
           legend=status.levels,
           pch=pch.status,
           title="Event type",
           bty='n')
    ## }}}
    
    ## {{{ Add warning
    if(add.warning){
        text(x=StudyLength/4,
             y=0.9*N,
             "This plot contains individual data. \nTherefore we cannot export it \noutside of Denmark Statistics!", pos=4, col="red", cex=1.5)
    }
    ## }}}

    ## {{{ inivsible output and seed
    .Random.seed <<- old
    out <- NULL
    ## browser()
    invisible(d[order(as.Date(d[,name.dateEntry])),])
    ## }}}
}


##----------------------------------------------------------------------
### PlotTimeDepTreat.R ends here

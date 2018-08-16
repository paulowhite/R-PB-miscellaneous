### permutation-test-function.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: May 16 2018 (09:11) 
## Version: 
## Last-Updated: Aug 16 2018 (15:37) 
##           By: Paul Blanche
##     Update #: 110
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


##' .. content for \description{Fully non-parametric permutation test to compare the longitudinal measurements (i.e. trajectory) of two groups} (no empty lines) ..
##'
##' .. content for \details{We permute the labels of the group variable. We compare a chi-square type distance computed with the actual data to the values observed under the null, obtain via permutation. Monte Carlo permutation are performed. A p-value is returned, which relates to the question: "Does a difference exist between the distributions of the longitudinal measurement of the two groups?". Plots help understanding the main ideas of what is computed. The function also returns: the distribution of the test statistics over the permuted data, the computation time, an error value and confidence interval, to evalute how precise is the estimate of the p-value, given the number of Monte Carlo runs.} ..
##' @title 
##' @param data data in long format
##' @param NMC number of monte carlo to approximate the permutation distribution (i.e. number of random permutations)
##' @param VarY name of outcome variable
##' @param VarGroup name of group variable (two groups to compare)
##' @param VarID name of the id variable, to identify obs from the same subject
##' @param varTime name of time variable
##' @param times the times that we want to keep (i.e. we do not look/compare to what is going on for other times )
##' @param showOne whether to visualize one random permuted dataset
##' @param seed to set the seed for reproducibility
##' @return 
##' @author Paul Blanche
##'
##' 

PermutTestCompare2Curves <- function(data,
                                     NMC=5000,
                                     VarY,
                                     VarGroup,
                                     VarID,
                                     varTime,
                                     times=NULL,
                                     showOne=TRUE,
                                     seed=2902){
    startkl <- Sys.time()
    ## browser()
    ## {{{ To plot p value with proper rounding
    specdec <- function(x, k) {
        xr <- round(x, k)
        if(xr==0){
            if(k>1){mystring <- paste0(paste(rep("0",k-1),collapse=""),1)
            }else{
                mystring <- "1"
            }
            paste0("<.",mystring)
        }else{
            format(xr, nsmall=k, scientific=FALSE)
        }
    }
    ## }}}
    ## {{{ Keep data with observation only at time points of interest
    if(!is.null(times)){
        data <- data[which(data[,varTime]%in% times),]
    }
    ## }}}
    ## browser()
    ## {{{ Test statistics realisation on actual data
    MeanGrT <- tapply(X=data[,VarY],INDEX=data[,c(VarGroup,varTime)],FUN=mean)
    VarGrT <- tapply(X=data[,VarY],INDEX=data[,c(VarGroup,varTime)],FUN=var)
    nGrT <- tapply(X=data[,VarY],INDEX=data[,c(VarGroup,varTime)],FUN=length)
    stat <- sum(((MeanGrT[1,]-MeanGrT[2,])^2)/(VarGrT[1,]/nGrT[1,] + VarGrT[2,]/nGrT[2,]))    
    ## }}}
    ## {{{ Prepare to permute
    set.seed(seed)
    allstats <- rep(NA,NMC)
    id <- unique(data[,VarID])   
    group <- sapply(id,function(x)data[which(data[,VarID]==x)[1],VarGroup])
    IDGr <- cbind(id,group)
    nValPersubj <- sapply(id,function(x)length(which(data[,VarID]==x)))
    ## }}}       
    ## {{{ permute and save test stat results
    for(i in 1:NMC){
        # shuffle data
        ds <- data
        NewGroup <- sample(group,size=length(group),replace=FALSE)
        ds[,VarGroup] <- rep(NewGroup,times=nValPersubj)
        MeanGrTs <- tapply(X=ds[,VarY],INDEX=ds[,c(VarGroup,varTime)],FUN=mean)
        VarGrTs <- tapply(X=ds[,VarY],INDEX=ds[,c(VarGroup,varTime)],FUN=var)
        nGrTs <- tapply(X=ds[,VarY],INDEX=ds[,c(VarGroup,varTime)],FUN=length)
        # test statistic
        allstats[i] <- sum(((MeanGrTs[1,]-MeanGrTs[2,])^2)/(VarGrTs[1,]/nGrTs[1,] + VarGrTs[2,]/nGrTs[2,])) 
    }
    p <- mean(allstats>stat)
    ForCIp <- binom.test(table(factor(allstats>stat,levels=c(TRUE,FALSE))))
    ## }}}   
    ## {{{ plot
    gr1 <- which(data[,VarGroup]==unique(group)[1])
    if(showOne){par(mfrow=c(1,3))}else{par(mfrow=c(1,2))}
    interaction.plot(data[gr1,varTime],
                     data[gr1,VarID],
                     data[gr1,VarY],
                     xlab="time",
                     ylab="Outcome Y",
                     legend=F,col="blue",
                     ylim=c(min(data[,VarY]),max(data[,VarY])))
    interaction.plot(data[-gr1,varTime],
                     data[-gr1,VarID],
                     data[-gr1,VarY],
                     legend=F,
                     add=TRUE,
                     col="red")
    title("Actual data")
    lines(1:ncol(MeanGrT),MeanGrT[1,],type="b",col="red",pch=16,lwd=3,lty=2,cex=3)
    lines(1:ncol(MeanGrT),MeanGrT[2,],type="b",col="blue",pch=16,lwd=3,lty=2,cex=3)
    legend("topleft",legend=paste0("group=",unique(group),", n=",c(sum(group==unique(group)[1]),sum(group==unique(group)[2]))),
           col=c("blue","red"),lwd=1,bty="n")
    legend("topright",
           legend=paste0("Pointwise mean in group=",unique(group)),
           col=c("blue","red"),#lwd=3,lty=2,
           pch=16, cex=1,
           bty="n")
    legend("bottomright",
           legend=paste0("Test stat=",specdec(stat,3)),
           pch=NA, 
           bty="n")
    if(showOne){
        gr1s <- which(ds[,VarGroup]==unique(group)[1])
        interaction.plot(ds[gr1s,varTime],
                         ds[gr1s,VarID],
                         ds[gr1s,VarY],
                         xlab="time",
                         ylab="Outcome Y",
                         legend=F,col="blue",
                         ylim=c(min(data[,VarY]),max(data[,VarY])))
        interaction.plot(ds[-gr1s,varTime],
                         ds[-gr1s,VarID],
                         ds[-gr1s,VarY],
                         legend=F,
                         add=TRUE,
                         col="red")
        lines(1:ncol(MeanGrTs),MeanGrTs[1,],type="b",col="red",pch=16,lwd=3,lty=2,cex=3)
        lines(1:ncol(MeanGrTs),MeanGrTs[2,],type="b",col="blue",pch=16,lwd=3,lty=2,cex=3)
        title("One randomly permuted dataset")
        legend("bottomright",
               legend=paste0("Test stat=",specdec(allstats[NMC],3)),
               pch=NA, 
               bty="n")
    }    
    ## browser()
    hist(allstats,xlim=range(1.1*c(allstats,stat)),
         main=paste0("Distribution under the null\n(",NMC,
                     " replications, p-value=",
                     specdec(p,3)," [",
                     specdec(ForCIp$conf.int[1],3),";",
                     specdec(ForCIp$conf.int[2],3),
                     "] )"),
         freq=FALSE,
         xlab="Test statistics")
    ## browser()
    abline(v=stat,col="red",lty=2)
    text(x=stat, y=0,pos=2,labels="Obs. \n test stat\n on actual data",col="red")
    ## }}}
    ## {{{ output
    stopkl <- Sys.time()    
    out <- list(pval=p,
                All=allstats,
                Stat=stat,
                CompTime=difftime(stopkl,startkl),
                NMC=NMC,
                error=1.96*sqrt(p*(1-p)/NMC),
                CIp=ForCIp$conf.int)
    ## }}}
    out
}

##----------------------------------------------------------------------
### permutation-test-function.R ends here

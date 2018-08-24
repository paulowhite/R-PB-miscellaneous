### permutation-test-function.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: May 16 2018 (09:11) 
## Version: 
## Last-Updated: Aug 24 2018 (15:35) 
##           By: Paul Blanche
##     Update #: 215
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
##' @param NMC number of monte carlo to approximate the permutation distribution (i.e. number of random permutations, if eaxct=FALSE)
##' @param VarY name of outcome variable
##' @param VarGroup name of group variable (two groups to compare)
##' @param VarID name of the id variable, to identify obs from the same subject
##' @param varTime name of time variable
##' @param times the times that we want to keep (i.e. we do not look/compare to what is going on for other times )
##' @param showOne whether to visualize one random permuted dataset
##' @param exact if TRUE, then all possible permutations are used (hence no Monte-Carlo approximation)
##' @param force if TRUE, it allows exact=TRUE even when this leads to more than 50.000 permutations. Otherwise it stops the computation and give a message error.
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
                                     exact=FALSE,
                                     force=FALSE,
                                     seed=2902){
    startkl <- Sys.time()
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
    ## {{{ Test statistics realisation on actual data
    MeanGrT <- tapply(X=data[,VarY],INDEX=data[,c(VarGroup,varTime)],FUN=mean)
    VarGrT <- tapply(X=data[,VarY],INDEX=data[,c(VarGroup,varTime)],FUN=var)
    nGrT <- tapply(X=data[,VarY],INDEX=data[,c(VarGroup,varTime)],FUN=length)
    stat <- sum(((MeanGrT[1,]-MeanGrT[2,])^2)/(VarGrT[1,]/nGrT[1,] + VarGrT[2,]/nGrT[2,]))    
    ## }}}
    ## {{{ Prepare to permute
    set.seed(seed)
    id <- unique(data[,VarID])   
    group <- sapply(id,function(x)data[which(data[,VarID]==x)[1],VarGroup])
    IDGr <- cbind(id,group)
    nValPersubj <- sapply(id,function(x)length(which(data[,VarID]==x)))
    ## {{{ Number of permutations: Monte Carlo number or number all possible  permutations if exact
    nSubj <- length(id)    
    if(exact){
        if(!force){
            NAllpermut <- factorial(nSubj)/(factorial(sum(group==unique(group)[1]))*factorial(sum(group==unique(group)[2])))
            if(NAllpermut>50000) stop("The use of exact=TRUE leads to computation for n=",NAllpermut," > 50,000 permuted data sets.\n Use force=TRUE if you really want to do such intensive computation.\n Otherwise, you should consider Monte-Carlo approximation using exact=FALSE and e.g. NMC=50000.\n")
        }
        AllCombPermut0 <- combn(x=nSubj,m=sum(group==unique(group)[1]))
        AllCombPermut <- apply(AllCombPermut0,2,
                               function(x){
                                   out <- rep(NA,nSubj)
                                   out[x] <- unique(group)[1]
                                   out[which(!(1:nSubj %in% x))] <- unique(group)[2]
                                   out
                               }
                               )
    }
    if(exact){
        nToLoop <- ncol(AllCombPermut)
    }else{
        nToLoop <- NMC
    }
    ## }}}
    allstats <- rep(NA,nToLoop)
    ## }}}
    ## {{{ permute and save test stat results
    for(i in 1:nToLoop){
        # shuffle data
        ds <- data
        if(!exact){
            NewGroup <- sample(group,size=length(group),replace=FALSE)
        }else{
            NewGroup <- AllCombPermut[,i]
        }
        ds[,VarGroup] <- rep(NewGroup,times=nValPersubj)
        MeanGrTs <- tapply(X=ds[,VarY],INDEX=ds[,c(VarGroup,varTime)],FUN=mean)
        VarGrTs <- tapply(X=ds[,VarY],INDEX=ds[,c(VarGroup,varTime)],FUN=var)
        nGrTs <- tapply(X=ds[,VarY],INDEX=ds[,c(VarGroup,varTime)],FUN=length)
        # test statistic
        allstats[i] <- sum(((MeanGrTs[1,]-MeanGrTs[2,])^2)/(VarGrTs[1,]/nGrTs[1,] + VarGrTs[2,]/nGrTs[2,]))
    }
    sumNA <- sum(is.na(allstats))
    if(sumNA>0){
        warning(paste0("N= ",sumNA," permuted data set out of ",nToLoop,
                       " (",specdec(sumNA/nToLoop*100,1),"%) lead to a test statistic to be NA. \n For those we set the value to be 1000 times the  value of the test statistic computed from the actual data, to ensure Type-I error control."))
        allstats[is.na(allstats)] <- 1000*stat
    }
    p <- mean(allstats>stat)
    if(!exact){ForCIp <- binom.test(table(factor(allstats>stat,levels=c(TRUE,FALSE))))}
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
    lines(1:ncol(MeanGrT),MeanGrT[unique(group)[1],],type="b",col="blue",pch=16,lwd=3,lty=2,cex=3)
    lines(1:ncol(MeanGrT),MeanGrT[unique(group)[2],],type="b",col="red",pch=16,lwd=3,lty=2,cex=3)
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
        lines(1:ncol(MeanGrTs),MeanGrTs[unique(group)[1],],type="b",col="blue",pch=16,lwd=3,lty=2,cex=3)
        lines(1:ncol(MeanGrTs),MeanGrTs[unique(group)[2],],type="b",col="red",pch=16,lwd=3,lty=2,cex=3)
        title("One randomly permuted dataset")
        legend("bottomright",
               legend=paste0("Test stat=",specdec(allstats[nToLoop],3)),
               pch=NA, 
               bty="n")
    }
    hist(allstats,xlim=range(1.1*c(allstats,stat)),
         main=paste0("Distribution under the null\n(",nToLoop,
                     " replications, p-value=",
                     specdec(p,3),
                     ifelse(exact,"",paste0(" [",
                                            specdec(ForCIp$conf.int[1],3),";",
                                            specdec(ForCIp$conf.int[2],3),
                                            "]")),
                     ")"),
         freq=FALSE,
         xlab="Test statistics")
    abline(v=stat,col="red",lty=2)
    text(x=stat, y=0,pos=2,labels="Obs. \n test stat\n on actual data",col="red")
    ## }}}
    ## {{{ output
    stopkl <- Sys.time()
    if(!exact){CIP <- ForCIp$conf.int}else{CIP <- NA}
    out <- list(pval=p,
                All=allstats,
                Stat=stat,
                CompTime=difftime(stopkl,startkl),
                NMC=NMC,
                nPermut=nToLoop,
                error=ifelse(!exact,1.96*sqrt(p*(1-p)/NMC),0),
                exact=exact,
                CIp=CIP)
    ## }}}
    out
}

##----------------------------------------------------------------------
### permutation-test-function.R ends here

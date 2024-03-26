#' Normal Simulation Function
#'
#' @param n sample size
#' @param p vector of percentiles of interest
#' @param stdev standard deviation of normal distribution
#' @param mean mean of normal distribution
#' @importFrom stats approxfun density dnorm integrate qnorm rnorm var
#'
#' @return absolute and relative bias of areas, means, variances, MSE, RSE
#' @export
#'
#' @examples \dontrun{normal.sim(3,c(.5,.6,.7,.8,.9,.95,.975),1,1)}
normal.sim <- function(n,p,stdev, mean){
  #setting up output
  k <- length(p)
  bias.area<-matrix(0,10000,k)
  area<-matrix(0,10000,k)
  means <- vector()
  vars <- vector()
  bias.mean<-vector()
  rel.bias.area <- vector()
  bias.var<-vector()
  removed <- 0
  MSE <- vector()
  SMSE <- vector()
  RSE <- vector()

  true.quantiles <- qnorm(p,mean,stdev) #quantiles based off of percentiles
  true.area <- 1-p #interested in area above quantile
  for(i in 1:10000){
    samp<-rnorm(n,mean,stdev)#draw new sample for each iteration

    #means
    means[i] <- mean(samp)
    bias.mean[i]<-abs(mean(samp)-mean)  #absolute bias
    avg.bias.mean <- mean(bias.mean)    #avg absolute bias
    rel.bias.mean <- avg.bias.mean/mean #relative bias of mean

    #variances
    vars[i] <- var(samp)
    bias.var[i]<-abs(var(samp)-stdev^2)    #var = std^2 absolute bias
    rel.bias.var <- mean(bias.var) / (stdev^2) #mean(abs bias)/ true var

    #integrating density estimate
    d<-density(samp, cut=5) #cut=5 should increase the upper and lower limits
    xx <- d$x
    yy <- d$y
    f <- approxfun(xx, yy,rule=2) #rule=2 says that for points outside interval, closest data extreme used
    p.unscaled<-vector()
    p.scaled<-vector()

    #allowing choice based on range of function
    lbound = min(mean-3*stdev,min(xx))
    ubound = max(mean+3*stdev,max(xx))

    C <- integrate(f,lbound,ubound, abs.tol=1e-40, subdivisions=100000, stop.on.error=FALSE)$value #entire density area
    for(j in 1:k){
      p.unscaled[j] <- integrate(f, true.quantiles[j], ubound, subdivisions=100000, abs.tol=1e-40,stop.on.error=FALSE)$value  #area above quantile
      p.scaled[j] <- p.unscaled[j] / C  #scale for proper probability

      if (p.scaled[j]>1){
        p.scaled[j] = NA #removing obviously incorrect calculations
        removed = removed + 1 #keeping track of number removed
      }
    }
    area[i,] <- p.scaled
    bias.area[i,]<-abs(p.scaled-true.area)

    #MSE
    true.d <- dnorm(xx,mean,stdev) #calculate heights of normal distribution at every point
    residuals <- yy - true.d #difference between estimate heights and actual heights
    MSE[i] <- mean(residuals^2)
    SMSE[i] <- MSE[i]/mean(yy^2) #standardized MSE = (obs-exp)^2/obs^2
    RSE[i] <- sum(residuals^2)/sum((mean(yy)-yy)^2) #RSE = (obs-exp)^2/(obs-mean)^2
  }

  #calculating summary stats
  avg.bias.area <- colMeans(bias.area,na.rm=TRUE)
  rel.bias.area <- avg.bias.area/true.area
  avg.var <- mean(vars)
  avg.bias.var <- mean(bias.var)
  T.stat <- avg.bias.mean/sqrt(avg.var/n)
  avg.MSE <- mean(MSE)
  avg.SMSE <- mean(SMSE)
  avg.RSE <- mean(RSE)


  list(avg.bias.area=avg.bias.area, rel.bias.area = rel.bias.area,
       avg.bias.mean=avg.bias.mean, rel.bias.mean=rel.bias.mean,
       rel.bias.var=rel.bias.var, avg.var=avg.var, avg.bias.var=avg.bias.var,
       T.stat=T.stat,
       removed=removed, avg.MSE=avg.MSE, avg.SMSE=avg.SMSE, avg.RSE=avg.RSE)
}

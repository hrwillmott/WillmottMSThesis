#' Gamma Simulation Function
#'
#' @param n sample size
#' @param p vector of percentiles of interest
#' @param alpha shape parameter for gamma function
#' @param beta scale parameter for gamma function
#' @importFrom stats approxfun density dgamma integrate qgamma rgamma var
#'
#' @return absolute and relative bias of areas, means, variances, MSE, RSE
#' @export
#'
#' @examples \dontrun{gamma.sim(3,c(.5,.6,.7,.8,.9,.95,.975),2,5)}
gamma.sim <- function(n,p,alpha,beta){
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

  true.quantiles <- qgamma(p,shape=alpha,scale=beta) #quantiles based off of percentiles
  true.area <- 1-p #interested in area above quantile
  for(i in 1:10000){
    samp<-rgamma(n,shape=alpha,scale=beta)#draw new sample for each iteration

    #means
    gam_mean <- alpha*beta #mean of gamma is alpha*beta
    means[i] <- mean(samp)
    bias.mean[i]<-abs(mean(samp)-gam_mean)  #absolute bias

    #variances
    gam_var <- alpha*beta^2 #var of gamma is alpha*beta^2
    vars[i] <- var(samp)
    bias.var[i]<-abs(var(samp)-gam_var)   #absolute bias

    #integrating density estimate
    d<-density(samp, cut=5) #cut=5 should increase the upper and lower limits
    xx <- d$x
    yy <- d$y
    f <- approxfun(xx, yy,rule=2) #rule=2 says that for points outside interval, closest data extreme used
    p.unscaled<-vector()
    p.scaled<-vector()

    #allowing choice based on range of function
    lbound = min(0,min(xx)) #gamma cannot go lower than 0
    ubound = max(gam_mean+3*sqrt(gam_var),max(xx)) #3std or max value of distribution

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
    posind <- which(xx>0) #only at values where density is evaluated above 0
    posx <- xx[posind]
    posy <- yy[posind]
    true.d <- dgamma(posx,shape=alpha,scale=beta) #calculate heights of gamma distribution at every point
    residuals <- posy - true.d #diff between est. heights and true heights
    MSE[i] <- mean(residuals^2)
    SMSE[i] <- MSE[i]/mean(posy^2) #standardized MSE = (obs-exp)^2/obs^2
    RSE[i] <- sum(residuals^2)/sum((mean(posy)-posy)^2) #RSE = (obs-exp)^2/(obs-mean)^2
  }

  #calculating summary stats
  avg.bias.area <- colMeans(bias.area,na.rm=TRUE)
  rel.bias.area <- avg.bias.area/true.area
  avg.var <- mean(vars)
  avg.bias.var <- mean(bias.var)
  rel.bias.var <- mean(bias.var) / (gam_var) #mean(abs bias)/ true var
  avg.bias.mean <- mean(bias.mean)    #avg absolute bias
  rel.bias.mean <- avg.bias.mean/gam_mean #relative bias of mean
  T.stat <- avg.bias.mean/sqrt(avg.var/n)
  avg.MSE <- mean(MSE)
  avg.SMSE <- mean(SMSE)
  avg.RSE <- mean(RSE)

  list(avg.bias.area=avg.bias.area, rel.bias.area = rel.bias.area,
       avg.bias.mean=avg.bias.mean, rel.bias.mean=rel.bias.mean,
       rel.bias.var=rel.bias.var, avg.var=avg.var,avg.bias.var=avg.bias.var,
       T.stat=T.stat,
       removed=removed, avg.MSE=avg.MSE, avg.SMSE=avg.SMSE, avg.RSE=avg.RSE)
}

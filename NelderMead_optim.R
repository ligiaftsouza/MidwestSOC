### Optimization using Nelder-Mead
# Arguments:
# d         : A data frame with columns "soc" and "depth"
# beta0     : An initial central guess for the beta decay parameter, various surround values also tried in separate optim call
# plotsv    : Default NA. Otherwise should be a filename (including path, no extension) where a plot is to be saved.

### Fits the function: soc = soc_f + (soc_i-soc_f)exp(-beta*depth) + epsilon
# epsilon is normal with mean zero, iid across depths
# soc_i is the value of the soc at depth 0
# soc_f is soc at max depth available.
# Actually transformed values of all parameters are fitted, to enforce the constraints
# soci, socf and beta should all be positive and socf should be less than soci.
# Multiple optimizations from different initial conditions surrounding these choices are run and the best result is returned.

nd_optim <-function(d,beta0,mkplot=FALSE)
{
  # The function to be optimized, depends on d, argument par has three parameters, call them x, y, z,
  # and they are unconstrained. Then soci=exp(x), socf=soci*(atan(y)+pi/2)/pi, beta=exp(z). 

  f<-function(par,d)
  {
    soci<-exp(par[1])
    socf<-soci*(atan(par[2])+pi/2)/pi
    beta<-exp(par[3])
    
    soc<-d$soc
    depth<-d$depth
    
    resids<-soc-socf-(soci-socf)*exp(-beta*depth)
    
    return(sum(resids^2))
  }
  
  # make sure the data are in depth order
  ord<-order(d$depth)
  d<-d[ord,]
  
  # now optimize
  rgseq<-seq(from=-3,to=3,length.out=13)
  # rgseq<-seq(from=-1,to=1,by=1)
  soci0<-d$soc[1]
  soci0<-ifelse(soci0<=0,1,soci0)
  socf0<-d$soc[dim(d)[1]]
  if (socf0>soci0)
  {
    socf0<-0.9*soci0
  }
  if (socf0==0)
  {
    socf0<-soci0/2
  }
  x0<-log(soci0)
  y0<-tan((socf0/soci0)*pi-pi/2)
  z0<-log(beta0)
  pars0<-as.matrix(expand.grid(x0+rgseq,y0+rgseq,z0+rgseq))

  wrap<-function(par)
  {
    res<-try(optim(par=par,fn=f,d=d,method="Nelder-Mead",
                   control=list(trace=0,maxit=5000)),silent=TRUE)
  }
  res<-apply(FUN=wrap,X=pars0,MARGIN=1)
  extractor<-function(x)
  {
    if (class(x)=="try-error")
    {
      return(NA)
    }else
    {
      return(x$value)
    }
  }
  allfs<-sapply(X=res,FUN=extractor)
  res<-res[[min(which(allfs==min(allfs)))]]
  
  # make the plot, if desired
  if (!is.na(mkplot) && class(res)!="try-error")
  {
    jpeg(file=paste0(mkplot,".jpg"))
    
    # the start parameters to fit
    depths<-seq(from=0,to=max(d$depth),length.out=100)
    midind<-ceiling(dim(pars0)[1]/2)
    pars<-pars0[midind,]
    soci<-exp(pars[1])
    socf<-soci*(atan(pars[2])+pi/2)/pi
    beta<-exp(pars[3])
    predsoc0<-socf+(soci-socf)*exp(-beta*depths)
    
    # the fit to plot
    pars<-res$par
    soci<-exp(pars[1])
    socf<-soci*(atan(pars[2])+pi/2)/pi
    beta<-exp(pars[3])
    predsoc<-socf+(soci-socf)*exp(-beta*depths)
    
    #now do the plot
    plot(d$depth,d$soc,type="o",xlab="depth",ylab="soc",ylim=range(d$soc,predsoc,predsoc0),
         main=paste("start=",f(pars0[midind,],d),"; end=",f(pars,d)))
    lines(depths,predsoc,type="l",col="green")
    lines(depths,predsoc0,type="l",col="red")
    
    dev.off()
  }
  
  # for later convenience, add the value of the objective function at the central set of start parems
  res$start_value<-f(c(x0,y0,z0),d)
  
  return(res)
}

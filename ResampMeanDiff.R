#Non-parametric, resampling-based alternative to a two-sample un-paired t-test.
#
#Args
#x        A set of values
#y        Another such
#num      The number of resamplings to do. Should be at least 1000, 10000 is better. 
#           Larger values give more accurate p-values.
#
#Output - a list with these named entries
#mnx      The mean of x
#mny      The mean of y
#p        The p-value for the test
#
#Details
#The difference mnx-mny is computed. The sets are resampled with replacement and 
#the same statistic is computed for each resampling. The p-value is the fraction of
#resampling-based values of the statistic which have the opposite sign as the 
#value computed on the original data.
#
ResampMeanDiff<-function(x,y,num)
{
  if (length(x)<=1 || length(y)<=1)
  {
    stop("Error in ResampMeanDiff: bad value for x or y")  
  }
  
  mnx<-mean(x)
  mny<-mean(y)
  mndiff<-mnx-mny
  
  mndiffre<-NA*numeric(num)
  for (counter in 1:num)
  {
    xre<-sample(x,replace=TRUE)
    yre<-sample(y,replace=TRUE)
    mndiffre[counter]<-mean(xre)-mean(yre)
  }
  
  p<-sum(sign(mndiffre)!=sign(mndiff))/num
  
  return(list(mnx=mnx,mny=mny,p=p))
}

# #Quick tests
# set.seed(101)
# x<-rnorm(100,0,1)
# y<-rnorm(100,.2,1)
# num<-10000
# ResampMeanDiff(x,y,num)
# 
# set.seed(101)
# x<-rnorm(100,0,1)
# y<-rnorm(100,.4,1)
# num<-10000
# ResampMeanDiff(x,y,num)

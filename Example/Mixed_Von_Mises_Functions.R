library('movMF')
library('date')
library('circular')
library('GoFKernel')
library('MCMCpack')
library('zoo')
library('LaplacesDemon')
den1=NULL
genden <- function(x,mu,kappa,com,K) {
  for (i in 1:K) {
    den1[i]=dvonmises(circular(x), mu=circular(mu[i]), kappa=kappa[i])*com[i]
  }
  return(sum(den1))
}
gendenv=Vectorize(genden,vectorize.args =c('x'))

den1=NULL
p <- function(x,mu,kappa,com,K) {
  for (i in 1:K) {
    den1[i]=pvonmises(circular(x), mu=circular(mu[i]), kappa=kappa[i],from = circular(0))*com[i]
  }
  return(sum(den1))
}
pv=Vectorize(p,vectorize.args =c('x'))

convert=function(c) {
  con1=as.numeric(unlist(strsplit(c, ":", fixed=TRUE)))
  con2=matrix(con1,nrow=length(con1)/3,ncol=3,byrow = TRUE)
  time=con2[,2]/60+con2[,1]
  return(time)
}
dbessel=function(x)
{
  if(x<709)
  {
    y=besselI(x, nu=1)/besselI(x, nu=0)
  } else
  {
    y=1
  }
  return(y)
}
cosdegree=function(cos,sin)
{
  if(sin<0)
    a=24-acos(cos)*24/(2*pi)
  else
    a=acos(cos)*24/(2*pi)
  return(a)
}
kesaifuc=list()
updatekesai=function(alpha,a,b,m,data,lambdabar)
{
  for (i in 1:n) {
    for (j in 1:K) {
      rou[i,j]=exp(digamma(alpha[j])-digamma(sum(alpha))-log(2*pi)
                   +a[j]/b[j]*t(m[j,])%*%data[i,]
                   -dbessel(lambdabar[j])*(a[j]/b[j]-lambdabar[j]))
    }
  }
  
  for (i in 1:n) {
    for (j in 1:K) {
      newkesai[i,j]=rou[i,j]/sum(rou[i,])
    }
  }
  if(sum(abs(newkesai-kesai))<0.001*n*K) 
  {kesaic=1}else 
    kesaic=0
  kesaifuc[[1]]=kesaic
  kesaifuc[[2]]=newkesai
  kesaifuc[[3]]=rou
  return(kesaifuc)
}
updatekesai2=function(alpha,a,b,m,data,lambdabar)
{
  for (i in 1:n) {
    for (j in 1:K) {
      rou[i,j]=exp(1/alpha[j]-1/sum(alpha)-log(2*pi)
                   +a[j]/b[j]*t(m[j,])%*%data[i,]
                   -besselI(lambdabar[j], nu=1)/besselI(lambdabar[j], nu=0)*(a[j]/b[j]-lambdabar[j]))
    }
  }
  
  for (i in 1:n) {
    for (j in 1:K) {
      newkesai[i,j]=rou[i,j]/sum(rou[i,])
    }
  }
  if(sum(abs(newkesai-kesai))<0.001*n*K) 
  {kesaic=1}else 
    kesaic=0
  kesaifuc[[1]]=kesaic
  kesaifuc[[2]]=newkesai
  return(kesaifuc)
}
alphafuc=list()
updatealpha=function(kesai)
{
  newalpha=alpha0+colSums(kesai)
  return(newalpha)
}
betamb=list()
newb=NULL
updatebetamb=function(kesai,data,lambdabar)
{
  newbeta=sqrt(diag((beta0*m0+t(kesai)%*%data)%*%t(beta0*m0+t(kesai)%*%data)))
  newm=(beta0*m0+t(kesai)%*%data)*(1/newbeta)
  
  for (i in 1:K) {
    newb[i]= b0[i]+besselI(lambdabar[i], nu=1)/besselI(lambdabar[i], nu=0)*sum(kesai[,i])+
      (dbessel(beta0[i]*lambdabar[i]))*beta0[i]
  }
  betamb[[1]]=newbeta
  betamb[[2]]=newm
  betamb[[3]]=newb
  return(betamb)
}
afuc=list()
newa=NULL
updatea=function(beta,lambdabar)
{
  for (i in 1:K) {
    newa[i]=a0[i]+beta[i]*lambdabar[i]*dbessel(beta[i]*lambdabar[i])
  }
  return(newa)
}
updatea2=function(beta,lambdabar)
{
  for (i in 1:K) {
    
    newa[i]=a0[i]+dbessel(beta[i]*lambdabar[i])
    
  }
  return(newa)
}
lambdabar=NULL
updatelambdabar=function(a,b)
{
  for (i in 1:K) {
    if(a[i]>1)
      lambdabar[i]=(a[i]-1)/b[i]
    else lambdabar[i]=a[i]/b[i]
  }
  return(lambdabar)
}









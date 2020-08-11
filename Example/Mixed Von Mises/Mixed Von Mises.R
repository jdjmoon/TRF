#read in the file 
s=read.csv('/Users/chentianyi/Downloads/sample.csv')
# save all time points into one vector for fitting
test=s[,c(4:7)]
time1=as.vector(t(test))
#get the number of time points
n=length(time1)
#convert time points into angle (0 to 2pi) then
#convert the angle into a 2d vector 
data=matrix(c(cos(time1*2*pi/24),sin(time1*2*pi/24)),ncol = 2, nrow = n,byrow = FALSE)
colnames(data)=c('COS','SIN')


## start fitting

# Algorithm is from 
#Taghia, Jalil & Ma, Zhanyu & Leijon, Arne. (2014). Bayesian Estimation of the von-Mises Fisher Mixture Model with Variational Inference. Pattern Analysis and Machine Intelligence, IEEE Transactions on. 36. 1701-1715. 10.1109/TPAMI.2014.2306426. 
#Object function here is to maximize lower bound.
#IN this paper authors don't have a way to decide the number of mixed components K and also the algorithm is subject to initialization difference.
#To solve these
#we used different K and chose the K with the max lower bound and different initializations to avoid the algorithm converge to local maximizer. 

resultlist=list()
LB=NULL  ## lower bound vector to record lower bound in each loop
kk=c(rep(2,3),rep(3,4),rep(4,6),rep(5,6))
for (f in 1:length(kk)) {
  K=kk[f]
  ##initial kesai 
  clusters <- kmeans(data[,1:2],K)
  table(clusters$cluster)
  kesai=matrix(0,nrow=n,ncol=K)
  for (i in 1:n) {
    kesai[i,clusters$cluster[i]]=1
  }
  #par(mfrow=c(1,2))
  #plot(data[,2],data[,1],col=clusters$cluster)
  #plot.circular(time,stack = T)
  ###initialize a,b,alpha,beta,lamda,m
  a=NULL
  b=NULL
  beta=NULL
  alpha=NULL
  lambdabar=NULL
  lb=m=NULL
  a0=rep(sample(c(.02,.01,.03),1),K)
  b0=rep(sample(c(.02,.01,.03),1),K)
  beta0=rep(sample(c(.02,.01,.03),1),K)
  alpha0=rep(sample(c(.5,.4,.6),1),K)
  lambdabar0=a0/b0
  m0=clusters$centers
  newkesai=rou=matrix(0,nrow = n,ncol = K)
  nn=20000
  c=10
  for (i in 1:nn) {
    if(i==1) {
      kesai=updatekesai(alpha0,a0,b0,m0,data,lambdabar0)[[2]]## all using initialization 
      alpha=updatealpha(kesai)
      beta=updatebetamb(kesai,data,lambdabar0)[[1]] ## 
      m=updatebetamb(kesai,data,lambdabar0)[[2]]
      b=updatebetamb(kesai,data,lambdabar0)[[3]]
      a=updatea(beta,lambdabar0)
      lambdabar=updatelambdabar(a,b)
      #check[i]=lambdabar
    } else
    {
      #print(c(updatekesai(alpha,a,b,m,data,lambdabar)[[1]],c))
      #print(updatekesai(alpha,a,b,m,data,lambdabar)[[3]][1,1])
      lb[i]=updatekesai(alpha,a,b,m,data,lambdabar)[[3]][1,1]
      if(updatekesai(alpha,a,b,m,data,lambdabar)[[1]]*c==1|lb[i]<10^(-5))
        break 
      kesai=updatekesai(alpha,a,b,m,data,lambdabar)[[2]]
      alpha=updatealpha(kesai)
      beta=updatebetamb(kesai,data,lambdabar)[[1]] ## 
      m=updatebetamb(kesai,data,lambdabar)[[2]]
      b=updatebetamb(kesai,data,lambdabar)[[3]]
      a=updatea(beta,lambdabar)
      if (sum(abs(updatelambdabar(a,b)-lambdabar))<0.001)
        c=1
      else c=0
      lambdabar=updatelambdabar(a,b)
    }
    #print(cosdegree(m[1,1],m[1,2]))
    #print(c)
  }
  
  LB[f]=lb[length(lb)]
  #print(LB[f])
  result=matrix(0,nrow = K,ncol=6)
  colnames(result)=c('Mean_Time','Concentration',"Proportion",'Beta','a','b')
  for (i in 1:K) {
    result[i,1]=cosdegree(m[i,1],m[i,2])
    result[i,2]=lambdabar[i]
    result[i,3]=alpha[i]
    result[i,4]=beta[i]
    result[i,5]=a[i]
    result[i,6]=b[i]
  }
  resultlist[[f]]=result
  print(f)
}
#get the max LB and corresponding fitting
best=which(LB==max(LB))[1]
result=resultlist[[best]]

##Figure 
#figure is range from 0 to 2*pi but is forced to be labeled as 0-24 for presentation purpose
hist(time1/24*2*pi,
     freq = F,col='green',
     ylim = c(0,1),
     xlim = c(0,2*pi),
     border = 'white',breaks = 50,xlab = 'Time',xaxt = "n",
     main = 'Original data vs fitting curve')
axis(side=1, at=seq(0,2*pi,length.out=25), labels=seq(0,24,by=1))
lines(x,gendenv(x,result[,1]/24*2*pi,result[,2],result[,3]/sum(result[,3]),nrow(result)),type='l',col='red')

##Error analysis

#MAE 
#using inverse method to sample 
#1, the inverse function of cumulative distribution function 
pvv=function(x) 
{
  y=pv(x,result[,1]/24*2*pi,result[,2],result[,3]/sum(result[,3]),nrow(result))
  return(y)
}
pinv=inverse(pvv,lower = 0,upper = 2*pi-0.000001)
pinvv=Vectorize(pinv)
s=runif(n)
sc=pinvv(s)
hist(sc,freq = FALSE,breaks = 30)
hist(time1/24*2*pi,breaks = seq(0,2*pi,length.out=48),freq = F)

#Calculate MAE
data=time1/24*2*pi
width=2*pi/50
breaks=seq(0,2*pi,width)
n=length(data)
MAE=NULL
for (i in 1:100) {
  s=runif(n)
  sc=pinvv(s)
  obcount=hist(data,breaks = breaks)$counts
  escounts=hist(sc,breaks = breaks)$counts
  MAE[i]=sum(abs(obcount-escounts))/length(escounts)
  print(i)
}
MAE
mean(MAE) #7.3504
quantile(MAE,0.025) #6.16 
quantile(MAE,0.975) #8.762 

#KLD
default=density(data,from = 0,to=2*pi)
estimated=gendenv(default$x,result[,1]/24*2*pi,result[,2],result[,3]/sum(result[,3]),nrow(result))
observed=default$y
plot(default$x,default$y,type = 'l',col='red ',ylim = c(0,1))
lines(default$x,estimated)
KLD(estimated,observed)$mean.sum.KLD




library(kernlab)
library(Rlab)
library(bdsmatrix)
library(MASS)
library(pracma)

set.seed(1)
x=runif(n=5000, min = -2*pi, max = 2*pi)
x
y=cos(x)+sin(x)
y


plot(x,y, xlab="time", ylab="Sin x + cos x")
n=length(x)
a=seq(-2*pi,2*pi,.1)
b=sin(a)+cos(a)
plot(a,b,type="l", xlab="time", ylab="Sin x + cos x")

########kernel

ker=function(n)
{
  start_r <- Sys.time()
  s=1
  l=1
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  
  I=diag(n)
  
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  
  
  inverse=solve(k+l*I)
  
  ff=y%*%inverse%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  plot(x,y)
  points(a,ff,col="blue")
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(cv,time_r)
  return(cv)
}


nn=c(50,100,250,500,750,1000,2000,3000,4000,5000)

try=ker(n=nn[1])

tt=matrix(NA,length(nn),3)
for(i in 1:length(nn)){
  ttt=ker(n=nn[i])
  tt[i,1]=nn[i]
  tt[i,2]=ttt[[1]]
  tt[i,3]=ttt[[2]]
  
}

tt1=tt
tt2=tt


tt1[10,3]=60*tt[10,3]
tt1


plot(nn,tt1[,3],xlab = "n",ylab="time in seconds")



#########RFF
rff=function(D,n)
{
  
  start_r <- Sys.time()
  s=1
  
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  y=cos(x)+sin(x)
  
  
  set.seed(1)
  w=rnorm(D,0,s)
  
  
  z2=matrix(NA,D,n)
  for (i in 1:D){
    for (j in 1:n){
      
      z2[i,j]=cos(w[i]*x[j])*(2/D)^(.5)
      
    }
  }
  
  
  z3=matrix(NA,D,n)
  for (i in 1:D){
    for (j in 1:n){
      
      z3[i,j]=sin(w[i]*x[j])*(2/D)^(.5)
      
    }
  }
  
  z4=rbind(z2,z3)
  z4=rbind(z4,y)
  z4=data.frame(t(z4))
  try2=lm(y~.-1,data=z4)
  summary(try2)
  
  
  cv=sum(((try2$fitted.values-y)^2)/length(try2$fitted.values))
  
  plot(x,try2$fitted.values)
  points(x,y,col="blue")
  
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  com=t(z3)%*%(z3)
  
  need=(k-com)
  pppf=norm(need,type="f")
  ppp1=norm(need,type="1")
  ppp2=norm(need,type="2")
  
  cv=list(cv,time_r,pppf,ppp1,ppp2)
  return(cv)
  
}
try1=rff(D=100,n=1000)



rfftt=matrix(NA,length(nn),6)
for(i in 1:length(nn)){
  ttt=rff(D=5,n=nn[i])
  rfftt[i,1]=nn[i]
  rfftt[i,2]=ttt[[1]]
  rfftt[i,3]=ttt[[2]]
  rfftt[i,4]=ttt[[3]]
  rfftt[i,5]=ttt[[4]]
  rfftt[i,6]=ttt[[5]]
  
  
}

rfftt1=rfftt
rfftt2=rfftt



plot(nn,rfftt1[,3],xlab = "n",ylab="time in seconds")














##########Sketching

ske=function(n,m)
{
  
  start_r <- Sys.time()
  s=1
  l=1
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  n=length(x)
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  I=diag(n)
  
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  set.seed(1)
  S=rnorm(n*m,0,1)
  
  S=matrix(S,n,m)
  
  term=t(S)%*%(k+l*I)%*%S
  
  w=solve(term)%*%t(S)%*%y
  al=S%*%w
  
  ff=(t(al)%*%(k1))
  cv=sum(((ff-b)^2)/length(ff))
  plot(x,y)
  points(a,ff,col="blue")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(cv,time_r)
  
  return(cv)
}

try2=ske(n=500,m=100)
try2

nn=c(50,100,250,500,750,1000,2000,3000,4000,5000)

try=ske(n=nn[1],m=5)

skett=matrix(NA,length(nn),3)
for(i in 1:length(nn)){
  ttt=ske(n=nn[i],m=max(50,nn[i]/10))
  skett[i,1]=nn[i]
  skett[i,2]=ttt[[1]]
  skett[i,3]=ttt[[2]]
  
}

skett1=skett
skett2=skett



plot(nn,skett1[,3],xlab = "n",ylab="time in seconds")


################
################
################
##########random walk
skeran=function(n,m)
{
  
  start_r <- Sys.time()
  s=1
  l=1
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  n=length(x)
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  I=diag(n)
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  ber=rbern(n*m, .5)
  S1=matrix(ber,n,m)
  S1<- ifelse(S1<0.1,-1,S1)
  
  term=t(S1)%*%(k+l*I)%*%S1
  
  w=solve(term)%*%t(S1)%*%y
  al=S1%*%w
  
  ff=(t(al)%*%(k1))
  cv=sum(((ff-b)^2)/length(ff))
  plot(x,y)
  points(a,ff,col="blue")
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(cv,time_r)
  return(cv)
  
}
try2tp=skeran(n=1000,m=100)

nn=c(50,100,250,500,750,1000,2000,3000,4000,5000)

try=skeran(n=nn[1],m=5)

skerantt=matrix(NA,length(nn),3)
for(i in 1:length(nn)){
  ttt=skeran(n=nn[i],m=max(50,nn[i]/10))
  skerantt[i,1]=nn[i]
  skerantt[i,2]=ttt[[1]]
  skerantt[i,3]=ttt[[2]]
  
}

skerantt1=skerantt
skerantt2=skerantt



plot(nn,skerantt1[,3],xlab = "n",ylab="time in seconds")




######Cholesky

ch=function(n)
{
  
  start_r <- Sys.time()
  s=1
  l=1
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  n=length(x)
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  
  I=diag(n)
  
  
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  
  r=chol(k+l*I)
  
  ri=backsolve(r, x = diag(ncol(r)))
  inv=ri%*%t(ri)
  ff=y%*%inv%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  plot(x,y)
  points(a,ff,col="blue")
  
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  
  
  cv=list(cv,time_r)
  return(cv)
}


chtt=matrix(NA,length(nn),3)
for(i in 1:length(nn)){
  ttt=ch(n=nn[i])
  chtt[i,1]=nn[i]
  chtt[i,2]=ttt[[1]]
  chtt[i,3]=ttt[[2]]
  
}

chtt1=chtt
chtt2=chtt


chtt1[9:10,3]=60*chtt[9:10,3]
chtt1


plot(nn,chtt1[,3],xlab = "n",ylab="time in seconds")







##########nystrom

ny=function(n,l,o)
{
  
  start_r <- Sys.time()
  s=1
  
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  n=length(x)
  
  
  ny1=seq(1,n,1)
  lo=sample(ny1,o)
  
  xo=x[lo]
  xno=x[-lo]
  
  w=matrix(NA,o,o)
  for(i in 1:o){
    for(j in 1:o){
      
      w[i,j]=exp(-s*((xo[i]-xo[j])^2)/2)
    }
  } 
  
  
  c1=matrix(NA,n-o,o)
  for(i in 1:n-o){
    for(j in 1:o){
      
      c1[i,j]=exp((xno[i]-xo[j])^2/2)
    }
  } 
  
  c=rbind(w,c1)
  wi=ginv(w)
  
  
  
  
  I=diag(n)
  
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  
  qi1=l*diag(1,o,o)
  qi2=wi%*%t(c)%*%c
  qi3=qi1+qi2
  qi4=ginv(qi3)
  qi5=c%*%qi4
  qi6=qi5%*%wi
  qi7=qi6%*%t(c)
  qi8=diag(1,n,n)
  inverse1=(qi8-qi7)/l
  
  #inverse=solve(kp+l*diag(1,n,n))
  #inverse2=solve(k+l*diag(1,n,n))
  
  ff=y%*%inverse1%*%k1
  #  ff1=y%*%inverse1%*%k1
  #  ff2=y%*%inverse2%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  #  cv1=sum(((ff1-b)^2)/length(ff1))
  #  cv2=sum(((ff2-b)^2)/length(ff2))
  
  plot(x,y)
  points(a,ff,col="blue")
  #  points(a,ff1,col="red")
  #  points(a,ff2,col="green")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  
  kp=c%*%wi%*%t(c)
  
  
  need=(k-kp)
  pppf=norm(need,type="f")
  ppp1=norm(need,type="1")
  ppp2=norm(need,type="2")
  
  cv=list(cv,time_r,pppf,ppp1,ppp2)
  return(cv)
  
}

try4=ny(n=1000,l = 125,o=10)
try4



nytt=matrix(NA,length(nn),6)
for(i in 1:length(nn)){
  ttt=ny(n=nn[i],l=nn[i]*1.25/10,o=10)
  nytt[i,1]=nn[i]
  nytt[i,2]=ttt[[1]]
  nytt[i,3]=ttt[[2]]
  nytt[i,4]=ttt[[3]]
  nytt[i,5]=ttt[[4]]
  nytt[i,6]=ttt[[5]]
  
  
}

nytt1=nytt




plot(nn,nytt1[,3],xlab = "n",ylab="time in seconds")







###############divide

ker1=function(n,a,x,y)
{
  start_r <- Sys.time()
  s=1
  l=1
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  
  I=diag(n)
  
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  
  
  inverse=solve(k+l*I)
  
  ff=y%*%inverse%*%k1
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(ff)
  return(cv)
}


dker=function(n,m)
{
  
  start_r <- Sys.time()
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  
  
  p=n/m
  
  data=data.frame(x,y)
  we=0
  for(i in 1:m){
    pdata=data[((i-1)*p+1):(i*p),]
    datax=pdata[,1]
    datay=pdata[,2]
    we[i]=ker1(n=p,a=a,x=datax,y=datay)
    
    
    
  }
  
  result=matrix(unlist(we),n,m)
  ff=apply(result,1,mean)
  
  
  cv=sum((ff-b)^2)/length(ff)
  
  plot(x,y)
  points(a,ff,col="blue")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  cv=list(cv,time_r)
  return(cv)
  
  
}





try111=dker(n=1000,m=10)




dtt=matrix(NA,length(nn),3)
for(i in 1:length(nn)){
  ttt=dker(n=nn[i],m=nn[i]/10)
  dtt[i,1]=nn[i]
  dtt[i,2]=ttt[[1]]
  dtt[i,3]=ttt[[2]]
  
}

dtt1=dtt



plot(nn,dtt1[,3],xlab = "n",ylab="time in seconds")






##############when n is 5000

fker=function(n)
{
  start_r <- Sys.time()
  s=1
  l=1
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  
  I=diag(n)
  
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  
  
  inverse=solve(k+l*I)
  
  ff=y%*%inverse%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  plot(x,y)
  points(a,ff,col="blue")
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(ff)
  return(cv)
}


fig1=fker(n=5000)




frff=function(D,n)
{
  
  start_r <- Sys.time()
  s=1
  
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  y=cos(x)+sin(x)
  
  
  set.seed(1)
  w=rnorm(D,0,s)
  
  
  z2=matrix(NA,D,n)
  for (i in 1:D){
    for (j in 1:n){
      
      z2[i,j]=cos(w[i]*x[j])*(2/D)^(.5)
      
    }
  }
  
  
  z3=matrix(NA,D,n)
  for (i in 1:D){
    for (j in 1:n){
      
      z3[i,j]=sin(w[i]*x[j])*(2/D)^(.5)
      
    }
  }
  
  z4=rbind(z2,z3)
  z4=rbind(z4,y)
  z4=data.frame(t(z4))
  try2=lm(y~.-1,data=z4)
  summary(try2)
  
  
  cv=sum(((try2$fitted.values-y)^2)/length(try2$fitted.values))
  
  plot(x,try2$fitted.values)
  points(x,y,col="blue")
  
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  cv=list(try2$fitted.values)
  return(cv)
  
}
fig2=frff(D=5,n=5000)




fske=function(n,m)
{
  
  start_r <- Sys.time()
  s=1
  l=1
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  n=length(x)
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  I=diag(n)
  
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  set.seed(1)
  S=rnorm(n*m,0,1)
  
  S=matrix(S,n,m)
  
  term=t(S)%*%(k+l*I)%*%S
  
  w=solve(term)%*%t(S)%*%y
  al=S%*%w
  
  ff=(t(al)%*%(k1))
  cv=sum(((ff-b)^2)/length(ff))
  plot(x,y)
  points(a,ff,col="blue")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(ff)
  
  return(cv)
}

fig3=fske(n=5000,m=500)




fskeran=function(n,m)
{
  
  start_r <- Sys.time()
  s=1
  l=1
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  n=length(x)
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  I=diag(n)
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  ber=rbern(n*m, .5)
  S1=matrix(ber,n,m)
  S1<- ifelse(S1<0.1,-1,S1)
  
  term=t(S1)%*%(k+l*I)%*%S1
  
  w=solve(term)%*%t(S1)%*%y
  al=S1%*%w
  
  ff=(t(al)%*%(k1))
  cv=sum(((ff-b)^2)/length(ff))
  plot(x,y)
  points(a,ff,col="blue")
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(ff)
  return(cv)
  
}
fig4=fskeran(n=5000,m=500)



fch=function(n)
{
  
  start_r <- Sys.time()
  s=1
  l=1
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  n=length(x)
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  
  I=diag(n)
  
  
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  
  r=chol(k+l*I)
  
  ri=backsolve(r, x = diag(ncol(r)))
  inv=ri%*%t(ri)
  ff=y%*%inv%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  plot(x,y)
  points(a,ff,col="blue")
  
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  
  
  cv=list(ff)
  return(cv)
}

fig5=fch(n=5000)


fny=function(n,l,o)
{
  
  start_r <- Sys.time()
  s=1
  
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  n=length(x)
  
  
  ny1=seq(1,n,1)
  lo=sample(ny1,o)
  
  xo=x[lo]
  xno=x[-lo]
  
  w=matrix(NA,o,o)
  for(i in 1:o){
    for(j in 1:o){
      
      w[i,j]=exp(-s*((xo[i]-xo[j])^2)/2)
    }
  } 
  
  
  c1=matrix(NA,n-o,o)
  for(i in 1:n-o){
    for(j in 1:o){
      
      c1[i,j]=exp((xno[i]-xo[j])^2/2)
    }
  } 
  
  c=rbind(w,c1)
  wi=ginv(w)
  
  
  
  
  I=diag(n)
  
  
  lx1=length(a)
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      k1[i,j]=exp(-s*((x[i]-x1[j])^2)/2)
    }
  }  
  
  qi1=l*diag(1,o,o)
  qi2=wi%*%t(c)%*%c
  qi3=qi1+qi2
  qi4=ginv(qi3)
  qi5=c%*%qi4
  qi6=qi5%*%wi
  qi7=qi6%*%t(c)
  qi8=diag(1,n,n)
  inverse1=(qi8-qi7)/l
  
  #inverse=solve(kp+l*diag(1,n,n))
  #inverse2=solve(k+l*diag(1,n,n))
  
  ff=y%*%inverse1%*%k1
  #  ff1=y%*%inverse1%*%k1
  #  ff2=y%*%inverse2%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  #  cv1=sum(((ff1-b)^2)/length(ff1))
  #  cv2=sum(((ff2-b)^2)/length(ff2))
  
  plot(x,y)
  points(a,ff,col="blue")
  #  points(a,ff1,col="red")
  #  points(a,ff2,col="green")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  cv=list(ff)
  return(cv)
  
}

fig5=fny(n=5000,l = 125*5,o=10)


fdker=function(n,m)
{
  
  start_r <- Sys.time()
  set.seed(1)
  x=runif(n, min = -2*pi, max = 2*pi)
  a=x
  y=cos(x)+sin(x)
  b=y
  
  
  
  p=n/m
  
  data=data.frame(x,y)
  we=0
  for(i in 1:m){
    pdata=data[((i-1)*p+1):(i*p),]
    datax=pdata[,1]
    datay=pdata[,2]
    we[i]=ker1(n=p,a=a,x=datax,y=datay)
    
    
    
  }
  
  result=matrix(unlist(we),n,m)
  ff=apply(result,1,mean)
  
  
  cv=sum((ff-b)^2)/length(ff)
  
  plot(x,y)
  points(a,ff,col="blue")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  cv=list(ff)
  return(cv)
  
  
}



fig6=fdker(n=5000,m=50)


plot(x,y,type="p",col="black",main="Comparison of different methods for n=5000"
     ,xlab="x",ylab="cos(x)+sin(x)",ylim=c(-1.5,3.5))
points(x,fig1[[1]],col="blue")
points(x,fig2[[1]],col="red")
points(x,fig3[[1]],col="green")
#points(x,fig4[[1]],col="yellow")
points(x,fig5[[1]],col="pink")
points(x,fig6[[1]],col="purple")
legend("topright",bty="n",legend=c("kernel","RFF","SKE","CHOL","NY","AVG"),lty=c(2,1.5),col=c("black","blue","red","green","pink","purple"))


####
plot(nn,chtt1[,2],col="green",type="l",xlab = "n",ylab="MSPE"
     ,main="MSPE for different methods",ylim=c(0,.32))
lines(nn,rfftt1[,2],col="blue")
lines(nn,skett1[,2],col="red")
lines(nn,tt1[,2],col="black")
lines(nn,nytt1[,2],col="pink")
lines(nn,dtt1[,2],col="purple")
legend("topright",bty="n",legend=c("kernel","RFF","SKE","CHOL","NY","AVG"),lty=c(2,2),col=c("black","blue","red","green","pink","purple"))


tab=c(chtt1[5,2],rfftt1[5,2],skett1[5,2],tt1[5,2],nytt1[5,2],dtt1[5,2])





######time

plot(nn,chtt1[,3],col="green",type="l",xlab = "n",ylab="time in seconds"
     ,main="Computation time for different methods")
lines(nn,rfftt1[,3],col="blue")
lines(nn,skett1[,3],col="red")
lines(nn,tt1[,3],col="black")
lines(nn,nytt1[,3],col="pink")
lines(nn,dtt1[,3],col="purple")
legend("topleft",bty="n",legend=c("kernel","RFF","SKE","CHOL","NY","AVG"),lty=c(10,80),col=c("black","blue","red","green","pink","purple"))




####error

plot(nn,chtt1[,2],col="green",type="l",xlab = "n",ylab="MSPE"
     ,main="MSPE for different methods",ylim=c(0,.32))
lines(nn,rfftt1[,2],col="blue")
lines(nn,skett1[,2],col="red")
lines(nn,tt1[,2],col="black")
lines(nn,nytt1[,2],col="pink")
lines(nn,dtt1[,2],col="purple")
legend("topright",bty="n",legend=c("kernel","RFF","SKE","CHOL","NY","AVG"),lty=c(2,2),col=c("black","blue","red","green","pink","purple"))




library(kernlab)
library(Rlab)
library(bdsmatrix)
library(MASS)
library(mvtnorm)
library(pracma)
library(coda)
library(data.table)
library(mcmcse)
library(truncnorm)
library(numDeriv)
library(MASS)
library(calibrate)
library(zipfR)





n=100
d=5
u=rep(0,d)
sigma=diag(1,d,d)
x=mvrnorm(n,u,sigma)

plot(x[,1],dmvnorm(x, u, sigma, log = FALSE))


ker=function(n)
{
  start_r <- Sys.time()
  
  d=5
  u=rep(0,d)
  sigma=diag(1,d,d)
  set.seed(1)
  x=mvrnorm(n,u,sigma)
  s=1
  l=2^(-5)
  a=x
  y=dmvnorm(x, u, sigma, log = FALSE)
  b=y
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x[j,h])^2/2)
      }
      k[i,j]=exp(-s*sum)
    }
  }  
  
  
  I=diag(n)
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  }  
  
  
  inverse=solve(k+l*I)
  
  ff=y%*%inverse%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(cv,time_r)
  return(cv)
}



try=ker(n=500)
try



nn=c(50,100,250,500,750,1000,1250,1500,1750,2000)



tt=matrix(NA,length(nn),3)
for(i in 1:length(nn)){
  ttt=ker(n=nn[i])
  tt[i,1]=nn[i]
  tt[i,2]=ttt[[1]]
  tt[i,3]=ttt[[2]]
  
}

tt1=tt
tt2=tt


tt1[9:10,3]=60*tt[9:10,3]
tt1


plot(nn,tt1[,3],xlab = "n",ylab="time in seconds")



##############RFF



rff=function(D,n,s)
{
  
  start_r <- Sys.time()
  
  
  set.seed(1)
  
  d=5
  u=rep(0,d)
  sigma=diag(1,d,d)
  set.seed(1)
  x=mvrnorm(n,u,sigma)
  
  a=x
  y=dmvnorm(x, u, sigma, log = FALSE)
  b=y
  
  
  set.seed(1)
  w=mvrnorm(D,u,s*sigma)
  
  
  z2=matrix(NA,D,n)
  for (i in 1:D){
    for (j in 1:n){
      
      z2[i,j]=sum(cos(w[i,]*x[j,]))*(2/D)^(.5)
      
    }
  }
  
  
  z3=matrix(NA,D,n)
  for (i in 1:D){
    for (j in 1:n){
      
      z3[i,j]=sum(sin(w[i,]*x[j,]))*(2/D)^(.5)
      
    }
  }
  
  z4=rbind(z2,z3)
  z4=rbind(z4,y)
  z4=data.frame(t(z4))
  try2=lm(y~.-1,data=z4)
  summary(try2)
  
  
  cv=sum((((try2$fitted.values)-y)^2)/length(try2$fitted.values))
  
  plot(x[,1],y)
  points(x[,1],(try2$fitted.values),col="blue")
  
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x[j,h])^2/2)
      }
      k[i,j]=exp(-s*sum)
    }
  }  
  
  com=t(z3)%*%(z3)
  
  need=(k-com)
  pppf=norm(need,type="f")
  ppp1=norm(need,type="1")
  ppp2=norm(need,type="2")
  
  cv=list(cv,time_r,pppf,ppp1,ppp2)
  return(cv)
  
}
try1=rff(D=150,n=300,s=175)
try1[[1]]


nn=c(50,100,250,500,750,1000,1250,1500,1750,2000)
dd=c(25,50,125,250,375,500,650,750,875,1000)
ss=c(1,10,100,450,1000,2000,3000,4500,6500,8500)



rfftt=matrix(NA,length(nn),6)
for(i in 1:length(nn)){
  
  ttt=rff(D=dd[i],n=nn[i],s=ss[i])
  rfftt[i,1]=nn[i]
  rfftt[i,2]=ttt[[1]]
  rfftt[i,3]=ttt[[2]]
  rfftt[i,4]=ttt[[3]]
  rfftt[i,5]=ttt[[4]]
  rfftt[i,6]=ttt[[5]]
  
  
}

rfftt1=rfftt
rfftt2=rfftt


plot(nn,rfftt1[,3],xlab = "n",ylab="time in seconds")


##############ske


ske=function(n,m)
{
  
  start_r <- Sys.time()
  
  d=5
  u=rep(0,d)
  sigma=diag(1,d,d)
  set.seed(1)
  x=mvrnorm(n,u,sigma)
  s=1
  l=2^(-5)
  a=x
  y=dmvnorm(x, u, sigma, log = FALSE)
  b=y
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x[j,h])^2/2)
      }
      k[i,j]=exp(-s*sum)
    }
  }    
  
  I=diag(n)
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  }  
  set.seed(1)
  S=rnorm(n*m,0,1)
  
  S=matrix(S,n,m)
  
  term=t(S)%*%(k+l*I)%*%S
  
  w=solve(term)%*%t(S)%*%y
  al=S%*%w
  
  ff=(t(al)%*%(k1))
  cv=sum(((ff-b)^2)/length(ff))
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(cv,time_r)
  
  return(cv)
}

try2=ske(n=500,m=300)
try2

skett=matrix(NA,length(nn),3)
for(i in 1:length(nn)){
  ttt=ske(n=nn[i],m=max(50,nn[i]/2))
  skett[i,1]=nn[i]
  skett[i,2]=ttt[[1]]
  skett[i,3]=ttt[[2]]
  
}

skett1=skett
skett2=skett
skett1[9:10,3]=60*skett[9:10,3]
skett1


plot(nn,skett1[,3],xlab = "n",ylab="time in seconds")





############chol



ch=function(n)
{
  start_r <- Sys.time()
  
  d=5
  u=rep(0,d)
  sigma=diag(1,d,d)
  set.seed(1)
  x=mvrnorm(n,u,sigma)
  s=1
  l=2^(-5)
  a=x
  y=dmvnorm(x, u, sigma, log = FALSE)
  b=y
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x[j,h])^2/2)
      }
      k[i,j]=exp(-s*sum)
    }
  } 
  
  I=diag(n)
  
  
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  } 
  
  r=chol(k+l*I)
  
  ri=backsolve(r, x = diag(ncol(r)))
  inv=ri%*%t(ri)
  ff=y%*%inv%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(cv,time_r)
  return(cv)
}

try3=ch(n=500)
try3





chtt=matrix(NA,length(nn),3)
for(i in 1:length(nn)){
  ttt=ch(n=nn[i])
  chtt[i,1]=nn[i]
  chtt[i,2]=ttt[[1]]
  chtt[i,3]=ttt[[2]]
  
}

chtt1=chtt
chtt2=chtt


chtt1[9:10,3]=60*chtt[9:10,3]
chtt1


plot(nn,chtt1[,3],xlab = "n",ylab="time in seconds")










############nystrom


ny=function(n,l,o)
{
  
  
  start_r <- Sys.time()
  
  d=5
  u=rep(0,d)
  sigma=diag(1,d,d)
  set.seed(1)
  x=mvrnorm(n,u,sigma)
  s=1
  
  a=x
  y=dmvnorm(x, u, sigma, log = FALSE)
  b=y
  
  
  I=diag(n)
  
  
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  } 
  
  
  ny1=seq(1,n,1)
  lo=sample(ny1,o)
  
  xo=x[lo,]
  xno=x[-lo,]
  
  
  w=matrix(NA,o,o)
  for(i in 1:o){
    for(j in 1:o){
      sum=0
      for(h in 1:5){
        sum= sum+((xo[i,h]-xo[j,h])^2/2)
      }
      w[i,j]=exp(-s*sum)
    }
  } 
  
  
  c1=matrix(NA,n-o,o)
  for(i in 1:n-o){
    for(j in 1:o){
      sum=0
      for(h in 1:5){
        sum= sum+((xno[i,h]-xo[j,h])^2/2)
      }
      c1[i,j]=exp(-s*sum)
    }
  } 
  
  
  
  c=rbind(w,c1)
  wi=ginv(w)
  kp=c%*%wi%*%t(c)
  
  
  qi1=l*diag(1,o,o)
  qi2=wi%*%t(c)%*%c
  qi3=qi1+qi2
  qi4=solve(qi3)
  qi5=c%*%qi4
  qi6=qi5%*%wi
  qi7=qi6%*%t(c)
  qi8=diag(1,n,n)
  inverse1=(qi8-qi7)/l
  
  #inverse=solve(kp+l*diag(1,n,n))
  #inverse2=solve(k+l*diag(1,n,n))
  
  ff=y%*%inverse1%*%k1
  #  ff1=y%*%inverse1%*%k1
  #  ff2=y%*%inverse2%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  #  cv1=sum(((ff1-b)^2)/length(ff1))
  #  cv2=sum(((ff2-b)^2)/length(ff2))
  
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  #  points(a,ff1,col="red")
  #  points(a,ff2,col="green")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  
  kp=c%*%wi%*%t(c)
  
  
  need=(k-kp)
  pppf=norm(need,type="f")
  ppp1=norm(need,type="1")
  ppp2=norm(need,type="2")
  
  cv=list(cv,time_r,pppf,ppp1,ppp2)
  return(cv)
  
}


nn=c(50,100,250,500,750,1000,1250,1500,1750,2000)
o=5
ll=c(5,8,17,37,53,77,103,120,140,140)

i=1
try4=ny(n=nn[i],l=ll[i],o=5)
try4




nytt=matrix(NA,length(nn),6)
for(i in 1:length(nn)){
  ttt=ny(n=nn[i],l=ll[i],o=5)
  nytt[i,1]=nn[i]
  nytt[i,2]=ttt[[1]]
  nytt[i,3]=ttt[[2]]
  nytt[i,4]=ttt[[3]]
  nytt[i,5]=ttt[[4]]
  nytt[i,6]=ttt[[5]]
  
  
}

nytt1=nytt




plot(nn,nytt1[,3],xlab = "n",ylab="time in seconds")


#######model avraging

ker1=function(n,x,a,y)
{
  start_r <- Sys.time()
  
  d=5
  u=rep(0,d)
  sigma=diag(1,d,d)
  
  s=1
  l=2^(-5)
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x[j,h])^2/2)
      }
      k[i,j]=exp(-s*sum)
    }
  }  
  
  
  I=diag(n)
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  }  
  
  
  inverse=solve(k+l*I)
  
  ff1=y%*%inverse
  ff=ff1%*%k1
  cv=ff
  return(cv)
}



dker=function(n,m)
{
  
  start_r <- Sys.time()
  
  d=5
  u=rep(0,d)
  sigma=diag(1,d,d)
  set.seed(1)
  x=mvrnorm(n,u,sigma)
  s=1
  l=2^(-5)
  a=x
  y=dmvnorm(x, u, sigma, log = FALSE)
  b=y
  data=data.frame(x,y)
  p=n/m
  
  data=data.frame(x,y)
  we=matrix(NA,n,m)
  for(i in 1:m){
    pdata=data[((i-1)*p+1):(i*p),]
    datax=pdata[,-6]
    datay=pdata[,6]
    we[,i]=ker1(n=p,x=datax,a=x,y=datay)
    
    
    
  }
  
  
  ff=apply(we,1,mean)
  
  
  cv=sum((ff-b)^2)/length(ff)
  
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  cv=list(cv,time_r)
  return(cv)
  
  
}


try111=dker(n=500,m=10)
try111


dtt=matrix(NA,length(nn),3)
for(i in 1:length(nn)){
  ttt=dker(n=nn[i],m=nn[i]/50)
  dtt[i,1]=nn[i]
  dtt[i,2]=ttt[[1]]
  dtt[i,3]=ttt[[2]]
  
}

dtt1[6:10,3]=60*dtt[6:10,3]
dtt1



plot(nn,dtt1[,3],xlab = "n",ylab="time in seconds")




######

######time

plot(nn,dtt1[,3],col="purple",type="l",xlab = "n",ylab="time in seconds"
     ,main="Computation time for different methods")
lines(nn,rfftt1[,3],col="blue")
lines(nn,skett1[,3],col="red")
lines(nn,tt1[,3],col="black")
lines(nn,nytt1[,3],col="pink")
lines(nn,chtt1[,3],col="green")
legend("topleft",bty="n",legend=c("kernel","RFF","SKE","CHOL","NY","AVG"),lty=c(10,80),col=c("black","blue","red","green","pink","purple"))




####error

plot(nn,nytt1[,2],col="pink",type="l",xlab = "n",ylab="MSPE"
     ,main="MSPE for different methods",ylim=c(0,3*10^(-07)))
lines(nn,rfftt1[,2],col="blue")
lines(nn,skett1[,2],col="red")
lines(nn,tt1[,2],col="black")
lines(nn,chtt1[,2],col="green")
lines(nn,dtt1[,2],col="purple")
legend("topright",bty="n",legend=c("kernel","RFF","SKE","CHOL","NY","AVG"),lty=c(2,2),col=c("black","blue","red","green","pink","purple"))














#########applied


adata=fread('https://archive.ics.uci.edu/ml/machine-learning-databases/parkinsons/telemonitoring/parkinsons_updrs.data')
head(adata)






y=adata[,6]
xx=adata[,-6]
x=xx[,-5]


a=x


aker=function(s,l)
{
  start_r <- Sys.time()
  
  y=adata[,6]
  y=data.frame(y)
  
  xx=adata[,-6]
  xx=xx[,-5]
  xx=xx[,-1]
  x=xx
  x=data.frame(x)
  a=x
  
  n=length(y[,1])
  b=y
  
  d=length(x[1,])
  
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      sum=0
      for(h in 1:d){
        sum= sum+((x[i,h]-x[j,h])^2/2)
      }
      k[i,j]=exp(-s*sum)
    }
  }  
  
  
  I=diag(n)
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:d){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  }  
  
  
  inverse=solve(k+l*I)
  
  ff=y%*%inverse%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(cv,time_r)
  return(cv)
}


try=aker(s=1,l=1)







arff=function(D,s)
{
  
  
  
  start_r <- Sys.time()
  
  y=adata[,6]
  
  
  xx=adata[,-6]
  xx=xx[,-5]
  xx=xx[,-1]
  x=xx
  x=data.frame(x)
  a=x
  
  n=length(y[,1])
  b=y
  
  d=length(x[1,])
  
  
  
  u=rep(0,d)
  sigma=diag(1,d,d)
  
  
  set.seed(1)
  w=mvrnorm(D,u,s*sigma)
  
  
  z2=matrix(NA,D,n)
  for (i in 1:D){
    for (j in 1:n){
      
      z2[i,j]=sum(cos(w[i,]*x[j,]))*(2/D)^(.5)
      
    }
  }
  
  
  z3=matrix(NA,D,n)
  for (i in 1:D){
    for (j in 1:n){
      
      z3[i,j]=sum(sin(w[i,]*x[j,]))*(2/D)^(.5)
      
    }
  }
  
  z4=rbind(z2,z3)
  z4=rbind(z4,y)
  z4=data.frame(t(z4))
  try2=lm(y~.-1,data=z4)
  summary(try2)
  
  
  cv=sum((((try2$fitted.values)-y)^2)/length(try2$fitted.values))
  
  plot(x[,1],y)
  points(x[,1],(try2$fitted.values),col="blue")
  
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  cv=list(cv,time_r)
  return(cv)
  
}

try1=arff(D=2500,s=1)




aske=function(m,s,l)
{
  
  start_r <- Sys.time()
  
  y=adata[,6]
  y=data.frame(y)
  
  xx=adata[,-6]
  xx=xx[,-5]
  xx=xx[,-1]
  x=xx
  x=data.frame(x)
  a=x
  
  n=length(y[,1])
  b=y
  
  d=length(x[1,])
  
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x[j,h])^2/2)
      }
      k[i,j]=exp(-s*sum)
    }
  }    
  
  I=diag(n)
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  }  
  set.seed(1)
  S=rnorm(n*m,0,1)
  
  S=matrix(S,n,m)
  
  term=t(S)%*%(k+l*I)%*%S
  
  w=solve(term)%*%t(S)%*%y
  al=S%*%w
  
  ff=(t(al)%*%(k1))
  cv=sum(((ff-b)^2)/length(ff))
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(cv,time_r)
  
  return(cv)
}

try2=aske(m=2500,s=1,l=1)





ach=function(s,l)
{
  start_r <- Sys.time()
  
  y=adata[,6]
  y=data.frame(y)
  
  xx=adata[,-6]
  xx=xx[,-5]
  xx=xx[,-1]
  x=xx
  x=data.frame(x)
  a=x
  
  n=length(y[,1])
  b=y
  
  d=length(x[1,])
  
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x[j,h])^2/2)
      }
      k[i,j]=exp(-s*sum)
    }
  } 
  
  I=diag(n)
  
  
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  } 
  
  r=chol(k+l*I)
  
  ri=backsolve(r, x = diag(ncol(r)))
  inv=ri%*%t(ri)
  ff=y%*%inv%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  cv=list(cv,time_r)
  return(cv)
}

try3=ach(s=1,l=1)



any=function(s,l,o)
{
  
  
  start_r <- Sys.time()
  
  y=adata[,6]
  y=data.frame(y)
  
  xx=adata[,-6]
  xx=xx[,-5]
  xx=xx[,-1]
  x=xx
  x=data.frame(x)
  a=x
  
  n=length(y[,1])
  b=y
  
  d=length(x[1,])
  
  I=diag(n)
  
  
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  } 
  
  
  ny1=seq(1,n,1)
  lo=sample(ny1,o)
  
  xo=x[lo,]
  xno=x[-lo,]
  
  
  w=matrix(NA,o,o)
  for(i in 1:o){
    for(j in 1:o){
      sum=0
      for(h in 1:5){
        sum= sum+((xo[i,h]-xo[j,h])^2/2)
      }
      w[i,j]=exp(-s*sum)
    }
  } 
  
  
  c1=matrix(NA,n-o,o)
  for(i in 1:n-o){
    for(j in 1:o){
      sum=0
      for(h in 1:5){
        sum= sum+((xno[i,h]-xo[j,h])^2/2)
      }
      c1[i,j]=exp(-s*sum)
    }
  } 
  
  
  
  c=rbind(w,c1)
  wi=ginv(w)
  kp=c%*%wi%*%t(c)
  
  
  qi1=l*diag(1,o,o)
  qi2=wi%*%t(c)%*%c
  qi3=qi1+qi2
  qi4=solve(qi3)
  qi5=c%*%qi4
  qi6=qi5%*%wi
  qi7=qi6%*%t(c)
  qi8=diag(1,n,n)
  inverse1=(qi8-qi7)/l
  
  #inverse=solve(kp+l*diag(1,n,n))
  #inverse2=solve(k+l*diag(1,n,n))
  
  ff=y%*%inverse1%*%k1
  #  ff1=y%*%inverse1%*%k1
  #  ff2=y%*%inverse2%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  #  cv1=sum(((ff1-b)^2)/length(ff1))
  #  cv2=sum(((ff2-b)^2)/length(ff2))
  
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  #  points(a,ff1,col="red")
  #  points(a,ff2,col="green")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*((x[i]-x[j])^2)/2)
    }
  }  
  
  
  kp=c%*%wi%*%t(c)
  
  
  need=(k-kp)
  pppf=norm(need,type="f")
  ppp1=norm(need,type="1")
  ppp2=norm(need,type="2")
  
  cv=list(cv,time_r,pppf,ppp1,ppp2)
  return(cv)
  
}


try4=any(s=1,l=700,o=5)




##define s and l
aker1=function(n,x,a,y)
{
  start_r <- Sys.time()
  
  
  n=length(y[,1])
  b=y
  
  d=length(x[1,])
  
  I=diag(n)
  
  
  u=rep(0,d)
  sigma=diag(1,d,d)
  
  s=1
  l=1
  
  
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x[j,h])^2/2)
      }
      k[i,j]=exp(-s*sum)
    }
  }  
  
  
  I=diag(n)
  
  
  lx1=length(a[,1])
  x1=a
  k1=matrix(NA,n,lx1)
  for(i in 1:n){
    for(j in 1:lx1){
      sum=0
      for(h in 1:5){
        sum= sum+((x[i,h]-x1[j,h])^2/2)
      }
      k1[i,j]=exp(-s*sum)
    }
  }  
  
  
  inverse=solve(k+l*I)
  
  ff1=y%*%inverse
  
  ff=ff1%*%k1
  cv=ff
  return(cv)
}



adker=function(m)
{
  
  start_r <- Sys.time()
  
  
  
  y=adata[,6]
  y=data.frame(y)
  
  xx=adata[,-6]
  xx=xx[,-5]
  xx=xx[,-1]
  x=xx
  x=data.frame(x)
  a=x
  
  n=length(y[,1])
  b=y
  
  d=length(x[1,])
  
  I=diag(n)
  
  
  s=1
  l=1
  
  
  data=data.frame(x,y)
  p=n/m
  
  data=data.frame(x,y)
  we=matrix(NA,n,m)
  for(i in 1:m){
    pdata=data[((i-1)*p+1):(i*p),]
    datax=pdata[,-(d+1)]
    datay=pdata[,(d+1)]
    we[,i]=ker1(n=p,x=datax,a=x,y=datay)
    
    
    
  }
  
  
  ff=apply(we,1,mean)
  
  
  cv=sum((ff-b)^2)/length(ff)
  
  plot(x[,1],y)
  points(a[,1],ff,col="blue")
  
  end_r <- Sys.time()
  time_r=end_r-start_r
  
  
  cv=list(cv,time_r)
  return(cv)
  
  
}


try111=dker(m=25)

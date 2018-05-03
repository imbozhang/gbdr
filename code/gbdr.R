################################################################################
# R Code for Groupwise Bayesian Dimension Reduction                                        
# Last updated: 2015/09/29                                                
################################################################################

rm(list=ls())
library(mvtnorm)
library(fields)   
library(dr)


################################################################################
# Simulation Model                                                             
################################################################################

#-------------------------------------------------------------------------------
# Model: Design 1 in Reich+2011                                                      
#-------------------------------------------------------------------------------

gen.design1 <-function(n, ps)
{
   ngrp<-length(ps)
   p<-sum(ps)

   b1<-c(1, 1, 1, 0, 0, 0)
   b2<-c(1, 0, 0, 0, 1, 3)
   b.true<-list(b1, b2)
   
   X<-matrix(rnorm(n*p),nrow=n)

   xb<-NULL
   
   for(g in 1:ngrp) { 
      X.g<-matrix(rnorm(n * ps[g]), nrow=n)
      xb.g<-as.vector(X %*% b.true[[g]])
      xb<-cbind(xb, xb.g)
   }
   
   V1<-xb[,1]
   V2 <- xb[,2]

   y<-0.4*V1^2 + 3*sin(0.25*V2) + 0.2*rnorm(n) 

   ans<-list(X=X, y=y, b.true=b.true, xb=xb)
   return(ans)
}

#-------------------------------------------------------------------------------
# Model: Design 2 in Reich+2011                                                      
#-------------------------------------------------------------------------------

gen.design2 <-function(n, ps)
{
   ngrp<-length(ps)
   p<-sum(ps)

   b1<-c(1, 1, 1, 0, 0, 0)
   b2<-c(1, 0, 0, 0, 1, 3)
   b.true<-list(b1, b2)

   X<-matrix(rnorm(n*p),nrow=n)

   xb<-NULL

   for(g in 1:ngrp) { 
      X.g<-matrix(rnorm(n * ps[g]), nrow=n)
      xb.g<-as.vector(X %*% b.true[[g]])
      xb<-cbind(xb, xb.g)
   }

   V1<-xb[,1]
   V2 <- xb[,2]

   y<-3*sin(0.25*V1) + 3*sin(0.25*V2) + 0.2*rnorm(n) 

   ans<-list(X=X, y=y, b.true=b.true, xb=xb)
   return(ans)
}

#-------------------------------------------------------------------------------
# Model (21) in groupwise dimension reduction                                                    
#-------------------------------------------------------------------------------

gen.data1<-function(n, ps)
{
   ngrp<-length(ps)
   p<-sum(ps)

   b1<-c(1, -1, 0, rep(0, ps[1]-3))
   b2<-c(1, 1, 1,  rep(0, ps[2]-3))
   b3<-c(1, 0, 0,  rep(0, ps[3]-3))
   b4<-c(0, 1, 0,  rep(0, ps[4]-3))
   b.true<-list(b1, b2, b3, b4)

   X<-xb<-NULL
   for(g in 1:ngrp) { 
      X.g<-matrix(rnorm(n * ps[g]), nrow=n)
      xb.g<-as.vector(X.g %*% b.true[[g]])
      X<-cbind(X, X.g)
      xb<-cbind(xb, xb.g)
   }

   y<-xb[,1] + 0.1*(xb[,2]+2)^2 + exp(0.5*xb[,3]) + sin(0.2*pi*xb[,4]) + 0.5*rnorm(n) 

   ans<-list(X=X, y=y, b.true=b.true, xb=xb)
   return(ans)
}

#-------------------------------------------------------------------------------
# Model (22) in groupwise dimension reduction                                                    
#-------------------------------------------------------------------------------

gen.data2<-function(n, ps)
{
   ngrp<-length(ps)
   p<-sum(ps)

   b1<-c(1, 1, 1, rep(0, ps[1]-3))
   b2<-c(1, 3, 0,  rep(0, ps[2]-3))
   b.true<-list(b1, b2)

   X<-xb<-NULL
   for(g in 1:ngrp) { 
      X.g<-matrix(rnorm(n * ps[g]), nrow=n)
      xb.g<-as.vector(X.g %*% b.true[[g]])
      X<-cbind(X, X.g)
      xb<-cbind(xb, xb.g)
   }

   y<-xb[,1] + xb[,1]*xb[,2] + 0.5*rnorm(n) 

   ans<-list(X=X, y=y, b.true=b.true, xb=xb)
   return(ans)
}

################################################################################
# Criterion                                                          
################################################################################


#-------------------------------------------------------------------------------
# LZC distance                                                                           
#-------------------------------------------------------------------------------

matpower = function(a,alpha){
if(ncol(a)==1) return(a^alpha) else
     small <- .00000000001
p1<-nrow(a)
eva<-eigen(a)$values
eve<-eigen(a)$vectors
eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
index<-(1:p1)[eva>small]
evai<-eva
evai[index]<-(eva[index])^(alpha)
ai<-eve%*%diag(evai)%*%t(eve)
return(ai)
}


dis = function(v1,v2){
p1 <- v1%*%matpower(t(v1)%*%v1,-1)%*%t(v1)
p2 <- v2%*%matpower(t(v2)%*%v2,-1)%*%t(v2)
d <- sum((p1-p2)*(p1-p2))
      return(d)
}

# Gram-Schmidt orthonormalization
orthnorm<-function(X)
{
   X<-as.matrix(X)
   n<-nrow(X)
   p<-ncol(X)

   W<-NULL
   if(p > 1) {
      W<-cbind(W, X[,1])
      for(k in 2:p) {
         gw<-rep(0, n)
         for(i in 1:(k-1)) {
            gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
            gw<-gw + gki * W[,i]
         }
         W<-cbind(W, X[,k] - gw)
      }
   } else {
      W<-cbind(W, X[,1])
   }

   W<-apply(W, 2, norm)
   W
}

# vector correlation and trace correlation between two spaces 
eval.space<-function(A, B, orthnm=TRUE) 
{
   if(!is.matrix(A)) A<-as.matrix(A)
   if(!is.matrix(B)) B<-as.matrix(B)
   if(orthnm) { 
      A<-orthnorm(A)
      B<-orthnorm(B) 
   }

   mat<-t(B) %*% A %*% t(A) %*% B
   d<-eigen(mat)$values
   d<-(d+abs(d))/2
   q<-sqrt(prod(d))
   r<-sqrt(mean(d))
   ans<-list(q=q, r=r)
   return(ans)
}

# evaluation
eval.d<-function(b.true, b.est) 
{
   ngrp<-min(length(b.true), length(b.est))
   ang<-vcr<-NULL
   for(g in 1:ngrp) {
      ang<-c(ang, angles(b.true[[g]], b.est[[g]]))
      vcr<-c(vcr, eval.space(b.true[[g]], b.est[[g]])$q)
   }
   ans<-rbind(ang, vcr)
   return(ans)
}

#-------------------------------------------------------------------------------
# angle between two spaces
#-------------------------------------------------------------------------------
angles<-function(B1, B2)
{
   if(!is.matrix(B1)) B1<-as.matrix(B1)
   if(!is.matrix(B2)) B2<-as.matrix(B2)

   if(ncol(B1) >= ncol(B2)) {
      B<-B1; B.hat<-B2
   } else {
      B<-B2; B.hat<-B1
   }

   P_1<-B %*% solve(t(B) %*% B) %*% t(B)
   if(ncol(B.hat) == 1) {
      nume<-as.vector(t(B.hat) %*% P_1 %*% B.hat)
      deno<-as.vector(t(B.hat) %*% B.hat)
      ratio<-nume / deno
   } else {
      BtB<-t(B.hat) %*% B.hat
      ei<-eigen(BtB)
      BtB2<-ei$vectors %*% diag(1/sqrt(ei$values)) %*% t(ei$vectors)
      M<-BtB2 %*% t(B.hat) %*% P_1 %*% B.hat %*% BtB2
      ratio<-abs(eigen(M)$values[nrow(M)])
   }
   ans<-acos(sqrt(ratio))/pi * 180
   if(ans > 90) ans<-180 - ans
   return(ans)
}


################################################################################
# Key function for BDR                                                      
################################################################################

g2z<-function(z,cuts,nc){
 g<-0*z
 for(c in 1:nc){g<-g+ifelse(z>cuts[c],1,0)}
g} 

rtnorm<-function(n,mu,sigma,L,U){
   l.prob<-pnorm(L,mu,sigma)
   u.prob<-pnorm(U,mu,sigma)
qnorm(runif(n,l.prob,u.prob),mu,sigma)}

Bayes.DR<-function(y,X,slices=7,d=1,ANOVA=T,
                   in.prob.a=1,in.prob.b=1,
                   iterations=5000,burnin=1000,update=100){
   small<-0.01
   eps<-0.1
   nc<-slices
   base<-1:d
   n<-length(y)
   p<-ncol(X)

   mean.x<-apply(X,2,mean,na.rm=T)
   sd.x<-apply(X,2,sd,na.rm=T)
   for(j in 1:p){
     X[,j]<-(X[,j]-mean.x[j])/sd.x[j]
   }
   y<-(y-mean(y,na.rm=T))/sd(y,na.rm=T)

   sd.cuts<-3*sqrt(p)
   cuts<-sd.cuts*qnorm(seq(0,1,length=nc+1))
   z<-matrix(rnorm(n*d),n,d)
   sir<-dr(y~X, method="sir", nslices=5)
   beta<-as.matrix(sir$evectors[,1:d])

   z<-X%*%beta

   inout<-prec.beta<-matrix(1,p,d)
   if(d>1){for(l in 2:d){
     beta[base[l-1],l:d]<-inout[base[l-1],l:d]<-0  
   }}
   
   prec.beta <- cbind(c(rep(1,3),rep(1000000,3)),c(rep(1000000,3),rep(1,3)))
   
   in.prob<-rep(1,d)
   g<-g2z(z,cuts,nc)

   mu<-array(rnorm(nc^d),rep(nc,d))
   taue<-1/var(y)
   tauz<-rep(1,d)
   int<-mean(y)
   taus<-1
   tXX<-t(X)%*%X

   mu1<-mu2<-mu3<-mu4<-rep(0,nc)
   tau1<-tau2<-tau3<-tau4<-1
   as<-bs<-20

   keep.beta<-array(0,c(iterations,p,d))
   inout.mn<-rep(0,p)

   for(i in 1:iterations){

     #fill in the zeros:
     int<-rnorm(1,mean(y-mu[g]),1/sqrt(taue*n))
     r<-y-int

     #update z
     Xbeta<-X%*%beta
     sige<-1/sqrt(taue)
     sigz<-1/sqrt(tauz)


     #update z
     for(l in 1:d){
       ppp<-matrix(0,n,nc)
       for(j in 1:nc){
         ggg<-g;ggg[,l]<-j
         ppp[,j]<-(pnorm(cuts[j+1],Xbeta[,l],sigz[l])-
                   pnorm(cuts[j],Xbeta[,l],sigz[l]))*
                   dnorm(r,mu[ggg],sige)
       }
       for(j in 1:n){
         if(!is.na(sum(ppp[j,]))){
           if(sum(ppp[j,])>0){
             g[j,l]<-sample(1:nc,1,prob=ppp[j,])
       }}}    

       z[,l]<-rtnorm(n,Xbeta[,l],sigz[l],cuts[g[,l]],cuts[g[,l]+1])
     }

     #update mu
     V<-M<-array(0,rep(nc,d))
     if(d==1){for(j in 1:nc){
       V[j]<-taue*sum(g==j)+taus*taue
       M[j]<-taue*sum(r[g==j])
     }}
     if(d==2){for(j1 in 1:nc){for(j2 in 1:nc){
       these<-(g[,1]==j1 & g[,2]==j2)
       V[j1,j2]<-taue*sum(these)+taue*taus
       M[j1,j2]<-taue*sum(r[these]) + taue*taus*(mu1[j1]+mu2[j2])
     }}}
     if(d==3){for(j1 in 1:nc){for(j2 in 1:nc){for(j3 in 1:nc){
       these<-(g[,1]==j1 & g[,2]==j2 & g[,3]==j3)
       V[j1,j2,j3]<-taue*sum(these)+taue*taus
       M[j1,j2,j3]<-taue*sum(r[these]) + 
                    taue*taus*(mu1[j1]+mu2[j2]+mu3[j3])
     }}}}
     if(d==4){for(j1 in 1:nc){for(j2 in 1:nc){for(j3 in 1:nc){for(j4 in 1:nc){
       these<-(g[,1]==j1 & g[,2]==j2 & g[,3]==j3 & g[,4]==j4)
       V[j1,j2,j3,j4]<-taue*sum(these)+taue*taus
       M[j1,j2,j3,j4]<-taue*sum(r[these]) + 
                    taue*taus*(mu1[j1]+mu2[j2]+mu3[j3]+mu[j4])
     }}}}}
     mu<-M/V+array(rnorm(nc^d),rep(nc,d))/sqrt(V)

     #update variances
     SS1<-sum((r-mu[g])^2)+tau1*sum(mu1^2)+tau2*sum(mu2^2)+tau3*sum(mu3^2)+tau4*sum(mu4^2)
     SS2<-0
     if(d==1){SS2<-sum(mu^2)}
     if(d==2){for(j1 in 1:nc){for(j2 in 1:nc){
       SS2<-SS2+(mu[j1,j2]-mu1[j1]-mu2[j2])^2
     }}}
     if(d==3){for(j1 in 1:nc){for(j2 in 1:nc){for(j3 in 1:nc){
       SS2<-SS2+(mu[j1,j2,j3]-mu1[j1]-mu2[j2]-mu3[j3])^2
     }}}}
     if(d==4){for(j1 in 1:nc){for(j2 in 1:nc){for(j3 in 1:nc){for(j4 in 1:nc){
       SS2<-SS2+(mu[j1,j2,j3,j4]-mu1[j1]-mu2[j2]-mu3[j3]-mu4[j4])^2
     }}}}}
     df<-n+nc^d+ifelse(ANOVA,d*nc,0)
     taue<-rgamma(1,df/2+eps,SS1/2+taus*SS2/2+eps)
     taus<-rgamma(1,(nc^d)/2+eps,taue*SS2/2+eps)

     if(!ANOVA){mu1<-mu2<-mu3<-mu4<-rep(0,nc)}
     if(d==2 & ANOVA){
       for(j in 1:nc){
         V<-1/(nc*taus*taue+taue*tau1)
         M<-sum(taus*taue*(mu[j,]-mu2))
         mu1[j]<-rnorm(1,V*M,sqrt(V))
     
         V<-1/(nc*taue*taus+taue*tau2)
         M<-sum(taus*taue*(mu[,j]-mu1))
         mu2[j]<-rnorm(1,V*M,sqrt(V))
       }
       tau1<-rgamma(1,nc/2+eps,taue*sum(mu1^2)/2+eps)
       tau2<-rgamma(1,nc/2+eps,taue*sum(mu2^2)/2+eps)
     }
   
     if(d==3 & ANOVA){
       for(j in 1:nc){
         V<-1/(nc*nc*taus*taue+taue*tau1)
         M<-0
         for(j1 in 1:nc){for(j2 in 1:nc){
           M<-M+taus*taue*(mu[j,j1,j2]-mu2[j1]-mu3[j2])
         }}
         mu1[j]<-rnorm(1,V*M,sqrt(V))
     
         V<-1/(nc*nc*taue*taus+taue*tau2)
         M<-0
         for(j1 in 1:nc){for(j2 in 1:nc){
           M<-M+taus*taue*(mu[j1,j,j2]-mu1[j1]-mu3[j2])
         }}
         mu2[j]<-rnorm(1,V*M,sqrt(V))

         V<-1/(nc*nc*taue*taus+taue*tau3)
         M<-0
         for(j1 in 1:nc){for(j2 in 1:nc){
           M<-M+taus*taue*(mu[j1,j2,j]-mu1[j1]-mu2[j2])
         }}
         mu3[j]<-rnorm(1,V*M,sqrt(V))
       }
       tau1<-rgamma(1,nc/2+eps,taue*sum(mu1^2)/2+eps)
       tau2<-rgamma(1,nc/2+eps,taue*sum(mu2^2)/2+eps)
       tau3<-rgamma(1,nc/2+eps,taue*sum(mu3^2)/2+eps)
       if(is.na(tau1)){tau1<-10}
       if(is.na(tau2)){tau2<-10}
       if(is.na(tau3)){tau3<-10}
     }
   
     if(d==4 & ANOVA){
       for(j in 1:nc){
         V<-1/(nc*nc*nc*taus*taue+taue*tau1)
         M<-0
         for(j1 in 1:nc){for(j2 in 1:nc){for(j3 in 1:nc){
           M<-M+taus*taue*(mu[j,j1,j2,j3]-mu2[j1]-mu3[j2]-mu4[j3])
         }}}
         mu1[j]<-rnorm(1,V*M,sqrt(V))
     
         V<-1/(nc*nc*nc*taue*taus+taue*tau2)
         M<-0
         for(j1 in 1:nc){for(j2 in 1:nc){for(j3 in 1:nc){
           M<-M+taus*taue*(mu[j1,j,j2,j3]-mu1[j1]-mu3[j2]-mu4[j3])
         }}}
         mu2[j]<-rnorm(1,V*M,sqrt(V))

         V<-1/(nc*nc*nc*taue*taus+taue*tau3)
         M<-0
         for(j1 in 1:nc){for(j2 in 1:nc){for(j3 in 1:nc){
           M<-M+taus*taue*(mu[j1,j2,j,j3]-mu1[j1]-mu2[j2]-mu4[j3])
         }}}
         mu3[j]<-rnorm(1,V*M,sqrt(V))

         V<-1/(nc*nc*nc*taue*taus+taue*tau3)
         M<-0
         for(j1 in 1:nc){for(j2 in 1:nc){for(j3 in 1:nc){
           M<-M+taus*taue*(mu[j1,j2,j,j3]-mu1[j1]-mu2[j2]-mu3[j3])
         }}}
         mu4[j]<-rnorm(1,V*M,sqrt(V))
       }
       tau1<-rgamma(1,nc/2+eps,taue*sum(mu1^2)/2+eps)
       tau2<-rgamma(1,nc/2+eps,taue*sum(mu2^2)/2+eps)
       tau3<-rgamma(1,nc/2+eps,taue*sum(mu3^2)/2+eps)
       tau4<-rgamma(1,nc/2+eps,taue*sum(mu4^2)/2+eps)
       if(is.na(tau1)){tau1<-10}
       if(is.na(tau2)){tau2<-10}
       if(is.na(tau3)){tau3<-10}
       if(is.na(tau4)){tau4<-10}
     }

     #update beta /tauz
     for(l in 1:d){       
       V<-solve(tauz[l]*tXX + diag(prec.beta[,l]))
       M<-tauz[l]*V%*%t(X)%*%z[,l]    
       if(i>5){beta[,l]<-M+t(chol(V))%*%rnorm(p)}
       big<-0*n*100+eps
       tauz[l]<-rgamma(1,n/2+eps,
                       sum((z[,l]-X%*%beta[,l])^2)/2+eps)
       if(is.na(tauz[l])){tauz[l]<-10}
       p1<-in.prob[l]*dnorm(beta[,l],0,1)
       p0<-(1-in.prob[l])*dnorm(beta[,l],0,small)
       inout[,l]<-rbinom(p,1,p1/(p0+p1))
       if(d>1){for(j in 2:d){inout[base[j-1],j:d]<-0}}
     }

     suc<-in.prob.a;fail<-in.prob.b
     for(j in 1:d){
       suc<-suc+sum(inout[,j])
       fail<-fail+sum(1-inout[,j])-(j-1)
     }
     in.prob[1:d]<-rbeta(1,suc,fail)

     if(i>burnin){
       inout.mn<-inout.mn+apply(inout,1,max)/(iterations-burnin)
     }

     keep.beta[i,,]<-beta

     if(i%%update==0){
       par(mfrow=c(1,1))
       if(d==1){ plot(beta,main=i)}
       if(d==2){ 
          plot(beta,col=0,main=i)
          text(beta[,1],beta[,2],1:p)
       }
       if(d>2){
         image.plot(1:d,1:p,t(beta),
         xlab="Direction",
         ylab="Covariate",main=i)
       }
     }
   }

  proj<-array(0,c(iterations,p,p))
  for(i in burnin:iterations){
    X<-keep.beta[i,,]
    proj[i,,]<-X%*%ginv(t(X)%*%X)%*%t(X)
  }
  P<-apply(proj[burnin:iterations,,],2:3,mean)

  
list(A=eigen(P)$vectors[,1:d],
     P=P,
     A.samples=keep.beta,
     in.prob=inout.mn
)}

################################################################################
# Simulation   											                                           
################################################################################

#-------------------------------------------------------------------------------
# Parameter Setting 
#-------------------------------------------------------------------------------
M=100                            # number of iterations of data sets
 
# data related
n <- 100
ps<-c(3, 3)
ds<-c(1, 1)
p1=ps[1]
p2=ps[2]
p=p1+p2
d1=ds[1]
d2=ds[2]
d=d1+d2

set.seed(100)
# result matrix
LZC = matrix(0,M,1)

#-------------------------------------------------------------------------------
# generate data from the example in BDR
#-------------------------------------------------------------------------------

ps<-c(3, 3)
ds<-c(1, 1)
p1=ps[1]
p2=ps[2]
p=p1+p2
d1=ds[1]
d2=ds[2]
d=d1+d2

for(h in 1:M){
data<-gen.design1(n=n, ps=ps)
X <- data$X
y <- data$y
b.true <- data$b.true



################################################################################
# BDR                                                             
################################################################################

out.BDR <- Bayes.DR(y,X, slices=7,d=d,ANOVA=T, in.prob.a=1,in.prob.b=1, iterations=5000,burnin=1000,update=100)
LZC[h] <- dis(cbind(b.true[[1]],b.true[[2]]), out.BDR$A)
}

result.mean <- mean(LZC)
result.sd <- sd(LZC)
result.mean
result.sd




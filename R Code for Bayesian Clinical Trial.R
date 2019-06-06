#The following is specified before enumeration algorithm.
sigma<-6
mu0<-4
mu1<-9
alpha<-0.1
beta<-0.1
gamma<-0.2
kappa0<-2
nu1<-mu0
nu2<-mu0
kappa1<-sigma/kappa0
kappa2<-sigma/kappa0
cT<-1-alpha
Nmax<-11
Nmin<-5
#Given the above specifications, we begin enumeration algorithm.
#n1 need to be set from n1 to n1+n2 in Hsiao(2008)
###############
#next, we should try different l1 & u1 values based on following value
l1domain<-c(2,3,4)
u1domain<-c(9,10,11)
step1<-matrix(NA,nrow=length(l1domain)*length(u1domain)*(Nmax-Nmin+1),ncol=4)
w<-1
for (i in Nmin:Nmax)
  for(j in 1:length(l1domain))
    for(k in 1:length(u1domain))
    {
      alphaf1<-1-pnorm(u1domain[k],mean=mu0,sd=sigma/sqrt(i))
      betaf1<-pnorm(l1domain[j],mean=mu1,sd=sigma/sqrt(i))
      gammaf<-1-pnorm(u1domain[k],mean=(mu0+mu1)/2,sd=sigma/sqrt(i))+pnorm(l1domain[j],mean=(mu0+mu1)/2,sd=sigma/sqrt(i))
      temp<-ifelse((alphaf1<alpha & betaf1<beta & gammaf<gamma),1,0)
      step1[w,]<-c(i,l1domain[j],u1domain[k],temp)
      w<-w+1
    }
newstep<-step1[step1[,4]==1,]
dim(newstep)
step2<-matrix(NA,nrow=nrow(newstep),ncol=5)
v<-1
n2<-31
  #next, we use monte carlo to compute alphaf2 and betaf2
  #first, compute alphaf2
  for (s in 1:nrow(newstep))
  {
   x1<-rnorm(n=10000,mean=mu0,sd=sigma/sqrt(newstep[s,1]))
   x2<-rnorm(n=10000,mean=mu0,sd=sigma/sqrt(n2))
   y2<-rnorm(n=10000,mean=mu0,sd=sigma/sqrt(n2))
   mu12<-sigma^2*nu1/(sigma^2+kappa1^2*(newstep[s,1]+n2))+(newstep[s,1]*x1+n2*x2)*kappa1^2/(sigma^2+kappa1^2*(newstep[s,1]+n2))
   tau12<-sqrt(1/((newstep[s,1]+n2)/sigma^2+1/kappa1^2))
   tau12<-rep(tau12,10000)
   mu2<-sigma^2*nu2/(sigma^2+kappa2^2*n2)+n2*y2*kappa2^2/(sigma^2+kappa2^2*n2)
   tau2<-sqrt(1/(n2/sigma^2+1/kappa2^2))
   tau2<-rep(tau2,10000)
   #compute PoP
   PoP<-rep(0,10000)
   for(i in 1:10000){
   integrand<-function(x) {dnorm(x,mean=mu12[i],sd=tau12[i])*pnorm(x,mean=mu2[i],sd=tau2[i])}
   PoP[i]<-integrate(integrand,lower=-Inf,upper=Inf)}
   PoP<-as.numeric(PoP)
   IndicatorofPoP<-rep(0,10000)
   for(i in 1:10000)
    {IndicatorofPoP[i]<-ifelse((PoP[i]>=cT),1,0)}
   Indicatorofx1<-rep(0,10000)
   for(i in 1:10000)
    {Indicatorofx1[i]<-ifelse((x1[i]>newstep[s,2] & x1[i]<newstep[s,3]),1,0)}
   alphaf2<-sum(IndicatorofPoP*Indicatorofx1)/10000
   #second, compute betaf2
   x1<-rnorm(n=10000,mean=mu1,sd=sigma/sqrt(newstep[s,1]))
   x2<-rnorm(n=10000,mean=mu1,sd=sigma/sqrt(n2))
   y2<-rnorm(n=10000,mean=mu0,sd=sigma/sqrt(n2))
   mu12<-sigma^2*nu1/(sigma^2+kappa1^2*(newstep[s,1]+n2))+(newstep[s,1]*x1+n2*x2)*kappa1^2/(sigma^2+kappa1^2*(newstep[s,1]+n2))
   tau12<-sqrt(1/((newstep[s,1]+n2)/sigma^2+1/kappa1^2))
   tau12<-rep(tau12,10000)
   mu2<-sigma^2*nu2/(sigma^2+kappa2^2*n2)+n2*y2*kappa2^2/(sigma^2+kappa2^2*n2)
   tau2<-sqrt(1/(n2/sigma^2+1/kappa2^2))
   tau2<-rep(tau2,10000)
   #compute PoP
   PoP<-rep(0,10000)
   for(i in 1:10000){
    integrand<-function(x) {dnorm(x,mean=mu12[i],sd=tau12[i])*pnorm(x,mean=mu2[i],sd=tau2[i])}
    PoP[i]<-integrate(integrand,lower=-Inf,upper=Inf)}
   PoP<-as.numeric(PoP)
   IndicatorofPoP<-rep(0,10000)
   for(i in 1:10000)
    {IndicatorofPoP[i]<-ifelse((PoP[i]<cT),1,0)}
   Indicatorofx1<-rep(0,10000)
   for(i in 1:10000)
    {Indicatorofx1[i]<-ifelse((x1[i]>newstep[s,2] & x1[i]<newstep[s,3]),1,0)}
   betaf2<-sum(IndicatorofPoP*Indicatorofx1)/10000
   alphaf<-alphaf1+alphaf2
   betaf<-betaf1+betaf2
   #next, we compute alphaB1, betaB1, alphastar1, betastar1
   x1<-rnorm(n=10000,mean=nu1,sd=sqrt(sigma^2/newstep[s,1]+kappa1^2))
   mu11<-sigma^2*nu1/(sigma^2+kappa1^2*newstep[s,1])+newstep[s,1]*x1*kappa1^2/(sigma^2+kappa1^2*newstep[s,1])
   tau11<-sqrt(1/(newstep[s,1]/sigma^2+1/kappa1^2))
   tau11<-rep(tau11,10000)
   part1<-pnorm(mu0,mean=mu11,sd=tau11)
   Indicatorofx1<-rep(0,10000)
   for(i in 1:10000)
    {Indicatorofx1[i]<-ifelse((x1[i]>newstep[s,3]),1,0)}
   alphaB1<-sum(Indicatorofx1*part1)/pnorm(mu0,mean=nu1,sd=kappa1)/10000
   newpart1<-1-pnorm(mu1,mean=mu11,sd=tau11)
   Indicatorofx1<-rep(0,10000)
   for(i in 1:10000)
    {Indicatorofx1[i]<-ifelse((x1[i]<newstep[s,2]),1,0)}
   betaB1<-sum(Indicatorofx1*newpart1)/(1-pnorm(mu1,mean=nu1,sd=kappa1))/10000
   alphastar1<-alphaB1*pnorm(mu0,mean=nu1,sd=kappa1)/(1-pnorm(newstep[s,3],mean=nu1,sd=sqrt(sigma^2/newstep[s,1]+kappa1^2)))
   betastar1<-betaB1*(1-pnorm(mu1,mean=nu1,sd=kappa1))/pnorm(newstep[s,2],mean=nu1,sd=sqrt(sigma^2/newstep[s,1]+kappa1^2))
   #next, we compute alphaB2, betaB2, alphastar2, betastar2
   x1<-rnorm(n=10000,mean=nu1,sd=sqrt(sigma^2/newstep[s,1]+kappa1^2))
   mu11<-sigma^2*nu1/(sigma^2+kappa1^2*newstep[s,1])+newstep[s,1]*x1*kappa1^2/(sigma^2+kappa1^2*newstep[s,1])
   tau11<-sqrt(1/(newstep[s,1]/sigma^2+1/kappa1^2))
   tau11<-rep(tau11,10000)
   x2<-rnorm(n=10000,mean=mu11,sd=sqrt(sigma^2/n2+tau11^2))
   y2<-rnorm(n=10000,mean=nu2,sd=sqrt(sigma^2/n2+kappa2^2))
   mu12<-sigma^2*nu1/(sigma^2+kappa1^2*(newstep[s,1]+n2))+(newstep[s,1]*x1+n2*x2)*kappa1^2/(sigma^2+kappa1^2*(newstep[s,1]+n2))
   tau12<-sqrt(1/((newstep[s,1]+n2)/sigma^2+1/kappa1^2))
   tau12<-rep(tau12,10000)
   mu2<-sigma^2*nu2/(sigma^2+kappa2^2*n2)+n2*y2*kappa2^2/(sigma^2+kappa2^2*n2)
   tau2<-sqrt(1/(n2/sigma^2+1/kappa2^2))
   tau2<-rep(tau2,10000)
   #compute PoP
   PoP<-rep(0,10000)
   for(i in 1:10000){
    integrand<-function(x) {dnorm(x,mean=mu12[i],sd=tau12[i])*pnorm(x,mean=mu2[i],sd=tau2[i])}
    PoP[i]<-integrate(integrand,lower=-Inf,upper=Inf)}
   PoP<-as.numeric(PoP)
   middlepart<-1-PoP
   IndicatorofPoP<-rep(0,10000)
   for(i in 1:10000)
    {IndicatorofPoP[i]<-ifelse((PoP[i]>=cT),1,0)}
   Indicatorofx1<-rep(0,10000)
   for(i in 1:10000)
    {Indicatorofx1[i]<-ifelse((x1[i]>newstep[s,2] & x1[i]<newstep[s,3]),1,0)}
   numeratorofalphaBandstar2<-sum(IndicatorofPoP*middlepart*Indicatorofx1)/10000
   denominatorofalphastar2<-sum(IndicatorofPoP*Indicatorofx1)/10000
   IndicatorofPoP<-rep(0,10000)
   for(i in 1:10000)
    {IndicatorofPoP[i]<-ifelse((PoP[i]<cT),1,0)}
   numeratorofbetaBandstar2<-sum(IndicatorofPoP*PoP*Indicatorofx1)/10000
   denominatorofbetastar2<-sum(IndicatorofPoP*Indicatorofx1)/10000
   #next, we derive the denominator of alphaB2 and betaB2
   x1<-rnorm(n=10000,mean=nu1,sd=sqrt(sigma^2/newstep[s,1]+kappa1^2))
   mu11<-sigma^2*nu1/(sigma^2+kappa1^2*newstep[s,1])+newstep[s,1]*x1*kappa1^2/(sigma^2+kappa1^2*newstep[s,1])
   tau11<-sqrt(1/(newstep[s,1]/sigma^2+1/kappa1^2))
   tau11<-rep(tau11,10000)
   #compute PofmuEmuS
   PofmuEmuS<-rep(0,10000)
   for(i in 1:10000){
    integrand<-function(x) {dnorm(x,mean=mu11[i],sd=tau11[i])*pnorm(x,mean=nu2,sd=kappa2)}
    PofmuEmuS[i]<-integrate(integrand,lower=-Inf,upper=Inf)}
   PofmuEmuS<-as.numeric(PofmuEmuS)
   Indicatorofx1<-rep(0,10000)
   for(i in 1:10000)
    {Indicatorofx1[i]<-ifelse((x1[i]>newstep[s,2] & x1[i]<newstep[s,3]),1,0)}
   denominatorofalphaB2<-sum((1-PofmuEmuS)*Indicatorofx1)/10000
   denominatorofbetaB2<-sum(PofmuEmuS*Indicatorofx1)/10000
   alphaB2<-numeratorofalphaBandstar2/denominatorofalphaB2
   betaB2<-numeratorofbetaBandstar2/denominatorofbetaB2
   alphastar2<-numeratorofalphaBandstar2/denominatorofalphastar2
   betastar2<-numeratorofbetaBandstar2/denominatorofbetastar2
   #next, we derive calibrator of alphaB2, betaB2,alphastar2, betastar2
   x1<-rnorm(n=10000,mean=nu1,sd=sqrt(sigma^2/newstep[s,1]+kappa1^2))
   mu11<-sigma^2*nu1/(sigma^2+kappa1^2*newstep[s,1])+newstep[s,1]*x1*kappa1^2/(sigma^2+kappa1^2*newstep[s,1])
   tau11<-sqrt(1/(newstep[s,1]/sigma^2+1/kappa1^2))
   tau11<-rep(tau11,10000)
   part1<-pnorm(mu0,mean=mu11,sd=tau11)
   Indicatorofx1<-rep(0,10000)
   for(i in 1:10000)
    {Indicatorofx1[i]<-ifelse((x1[i]>newstep[s,2] & x1[i]<newstep[s,3]),1,0)}
   calofalphaB2<-sum(Indicatorofx1*part1)/pnorm(mu0,mean=nu1,sd=kappa1)/10000
   newpart1<-1-pnorm(mu1,mean=mu11,sd=tau11)
   calofbetaB2<-sum(Indicatorofx1*newpart1)/(1-pnorm(mu1,mean=nu1,sd=kappa1))/10000
   newnewpart1<-pnorm(mu1,mean=mu11,sd=tau11)-pnorm(mu0,mean=mu11,sd=tau11)
   Indicatorofx1<-rep(0,10000)
   for(i in 1:10000)
    {Indicatorofx1[i]<-ifelse((x1[i]>newstep[s,3]),1,0)}
   calofalphastar2<-sum(newnewpart1*Indicatorofx1)/(1-pnorm(newstep[s,3],mean=nu1,sd=sqrt(sigma^2/newstep[s,1]+kappa1^2)))/10000
   Indicatorofx1<-rep(0,10000)
   for(i in 1:10000)
    {Indicatorofx1[i]<-ifelse((x1[i]<newstep[s,2]),1,0)}
   calofbetastar2<-sum(newnewpart1*Indicatorofx1)/(1-pnorm(newstep[s,2],mean=nu1,sd=sqrt(sigma^2/newstep[s,1]+kappa1^2)))/10000
   alphaB2cal<-alphaB2*calofalphaB2
   betaB2cal<-betaB2*calofbetaB2
   alphastar2cal<-alphastar2*calofalphastar2
   betastar2cal<-betastar2*calofbetastar2
   alphaB<-alphaB1+alphaB2cal
   betaB<-betaB1+betaB2cal
   alphastar<-alphastar1+alphastar2cal
   betastar<-betastar1+betastar2cal
   result<-ifelse((alphaf<alpha & betaf<beta & alphaB<alpha & betaB<beta & alphastar<alpha & betastar<beta),1,0)
   step2[v,]<-c(newstep[s,1:3],n2,result)
   v<-v+1
  }
newnewstep<-step2[step2[,5]==1,]
dim(newnewstep)
ESSB<-newnewstep[,1]+2*newnewstep[,4]*(pnorm(newnewstep[,3],mean=nu1,sd=sqrt(sigma^2/newnewstep[,1]+kappa1^2))-pnorm(newnewstep[,2],mean=nu1,sd=sqrt(sigma^2/newnewstep[,1]+kappa1^2)))
min(ESSB)
which.min(ESSB)

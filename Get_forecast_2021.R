#### CLEAR MEMORY
rm(list=ls())
#### Prevent Warnings
options(warn=-1)
#### Clear the R console
cat("\f")




# read data ---------------------------------------------------------------
Year=2021

time_elapsed=0
ptm1 = proc.time()
mainDir=getwd()
knitr::opts_knit$set(root.dir = mainDir)
setwd(mainDir)
suppressPackageStartupMessages(source("library.R"))
source("library.R")

load_packages(c("parallel","rstan","ggplot2","rlist","sf","raster","tiff","grid","MVN","knitr","FNN","foreach","doParallel","seqinr"),quietly=TRUE)
name_dir="Results"
dir_create("Results")
dir_r=paste(mainDir,"Results",sep = "/")
dir_create("Results/Parameters")
dir_data_r=paste(mainDir,"Results/Parameters",sep = "/")
dir_create("Results/Season2021")
dir_JPG=paste(dir_r,"Season2021",sep = "/")

setwd(mainDir)
Q_F=list.load(paste('Data_Covariates_1to5day_LeadTimes_',Year,'.rds',sep=""))
data=list.load("Observed_data_forecast_2021.rds")
Mod=list(desc="")
lt_ind=seq(30,26,-1)
lt_ind2=seq(30,26,-1)
for (j in 1:5) {
  #j-day lead
  setwd(dir_data_r)
  Mod[[paste(j,"_d",sep="")]]=list.load(paste("Par_bayesian_model_gamma_best_mod_",j,"LT.rds",sep=""))
  
}


N2=dim(Q_F$VIC_For_1d)[1]
N=dim(Mod$`1_d`$beta_Gar)[1]
tmp=array(NA,dim=c(N,N2,5),dimnames = list(1:N,1:N2,colnames(Q_F$Qthreshold)))
for (i in 1:5) {
  Q_F[[paste("Q_",i,"d",sep="")]]=tmp
}


# get forecast only if Model_fitted=1 -------------------------------------

##1 day lead
k=1
datath=Q_F[[paste("VIC_For_",k,"d",sep="")]][,4:8]
tmp1=Q_F[[paste("VIC_Sim_",k,"d",sep="")]][,5:9]
tmp2=Q_F[[paste("P_Ob_",k,"d",sep="")]][,4:8]
covP=as.matrix(cbind(tmp1[,1:2],tmp2[,3],tmp1[,4:5]))
Qth=Mod[[paste(k,"_d",sep="")]]$thrQ
mod_1=Mod[[paste(k,"_d",sep="")]]
#get the element for each step
ind_f=list(desc="")
nam=c("Gar","Man","Han","Hos","San")
for (i in 1:5) {
  ind_f[[paste("ind_",nam[i],1,sep="")]]=which(datath[,i]<=Qth[i])
  ind_f[[paste("N_",nam[i],1,sep="")]]=length(ind_f[[paste("ind_",nam[i],1,sep="")]])
  ind_f[[paste("ind_",nam[i],2,sep="")]]=which(datath[,i]>Qth[i])
  ind_f[[paste("N_",nam[i],2,sep="")]]=length(ind_f[[paste("ind_",nam[i],2,sep="")]])
  
}
exp1=c(2,2.5,1.5,2,2)

time_elapsed=0
ptm = proc.time()
numCores <- detectCores()-1
registerDoParallel(numCores)
for (i in 1:5) {
  res1=foreach (j=1:N2, .packages=c("stats")) %dopar% {
    beta=mod_1[[paste("beta_",nam[i],sep = "")]]
    phi=mod_1[[paste("phi_",nam[i],sep = "")]]
    sim=rep(NA,dim(beta)[1])
    ind_t=which(ind_f[[paste("ind_",nam[i],1,sep="")]]==j)
    if(length(ind_t)>0){
      mu=beta[,1]+beta[,2]*datath[j,i]^exp1[i]+beta[,3]*covP[j,i]
      sigma = phi[,1]+phi[,2]*datath[j,i]+phi[,3]*covP[j,i]
    }else{
      mu=beta[,4]+beta[,5]*datath[j,i]^exp1[i]+beta[,6]*covP[j,i]
      sigma = phi[,4]+phi[,5]*datath[j,i]+phi[,6]*covP[j,i]
    }
    # if(i==3){
    #   ind_t2=which(ind_f[[paste("ind_",nam[i],2,sep="")]]==j)
    #   ind_t3=which(ind_f[[paste("ind_",nam[i],3,sep="")]]==j)
    #   if(length(ind_t)>0){
    #     mu=beta[,1]+beta[,2]*datath[j,i]^exp1[i]
    #     sigma = phi[,1]+phi[,2]*datath[j,i]
    #   }else if (length(ind_t2)>0){
    #     mu=beta[,3]+beta[,4]*datath[j,i]^exp1[i]+beta[,5]*covP[j,i]
    #     sigma = phi[,3]+phi[,4]*datath[j,i]+phi[,5]*covP[j,i]
    #   }else if (length(ind_t3)>0){
    #     mu=beta[,6]+beta[,7]*datath[j,i]^exp1[i]+beta[,8]*covP[j,i]
    #     sigma = phi[,6]+phi[,7]*datath[j,i]+phi[,8]*covP[j,i]
    #   }
    # }else{
    #   if(length(ind_t)>0){
    #     mu=beta[,1]+beta[,2]*datath[j,i]^exp1[i]+beta[,3]*covP[j,i]
    #     sigma = phi[,1]+phi[,2]*datath[j,i]+phi[,3]*covP[j,i]
    #   }else{
    #     mu=beta[,4]+beta[,5]*datath[j,i]^exp1[i]+beta[,6]*covP[j,i]
    #     sigma = phi[,4]+phi[,5]*datath[j,i]+phi[,6]*covP[j,i]
    #   }
    # }
    
    
    shape = mu^2/sigma^2
    rate = mu/sigma^2
    for (k in 1:N) {
      sim[k]=rgamma(1,shape = shape[k],rate = rate[k])
    }
    sim
  }
  Q_F[[paste("Q_",k,"d",sep="")]][,,i]=t(matrix(unlist(res1), ncol = N, byrow = T))
  
}
stopImplicitCluster()


##2 day lead

k=2
datath=Q_F[[paste("VIC_For_",k,"d",sep="")]][,4:8]
tmp1=Q_F[[paste("VIC_Sim_",k,"d",sep="")]][,5:9]
covP=tmp1
Qth=Mod[[paste(k,"_d",sep="")]]$thrQ
mod_1=Mod[[paste(k,"_d",sep="")]]
#get the element for each step
ind_f=list(desc="")
nam=c("Gar","Man","Han","Hos","San")
for (i in 1:5) {
  ind_f[[paste("ind_",nam[i],1,sep="")]]=which(datath[,i]<=Qth[i])
  ind_f[[paste("N_",nam[i],1,sep="")]]=length(ind_f[[paste("ind_",nam[i],1,sep="")]])
  ind_f[[paste("ind_",nam[i],2,sep="")]]=which(datath[,i]>Qth[i])
  ind_f[[paste("N_",nam[i],2,sep="")]]=length(ind_f[[paste("ind_",nam[i],2,sep="")]])
  
}
exp1=c(1.5,2,1.5,1.5,1.5)

time_elapsed=0
ptm = proc.time()
numCores <- detectCores()-1
registerDoParallel(numCores)
for (i in 1:5) {
  res1=foreach (j=1:N2, .packages=c("stats")) %dopar% {
    beta=mod_1[[paste("beta_",nam[i],sep = "")]]
    phi=mod_1[[paste("phi_",nam[i],sep = "")]]
    sim=rep(NA,dim(beta)[1])
    ind_t=which(ind_f[[paste("ind_",nam[i],1,sep="")]]==j)
    if(length(ind_t)>0){
      mu=beta[,1]+beta[,2]*datath[j,i]^exp1[i]+beta[,3]*covP[j,i]
      sigma = phi[,1]+phi[,2]*datath[j,i]+phi[,3]*covP[j,i]
    }else{
      mu=beta[,4]+beta[,5]*datath[j,i]^exp1[i]+beta[,6]*covP[j,i]
      sigma = phi[,4]+phi[,5]*datath[j,i]+phi[,6]*covP[j,i]
    }
    
    shape = mu^2/sigma^2
    rate = mu/sigma^2
    for (k in 1:N) {
      sim[k]=rgamma(1,shape = shape[k],rate = rate[k])
    }
    sim
  }
  Q_F[[paste("Q_",k,"d",sep="")]][,,i]=t(matrix(unlist(res1), ncol = N, byrow = T))
  
}
stopImplicitCluster()



##3 day lead

k=3
datath=Q_F[[paste("VIC_For_",k,"d",sep="")]][,4:8]
tmp1=Q_F[[paste("VIC_Sim_",k,"d",sep="")]][,5:9]
covP=tmp1
Qth=Mod[[paste(k,"_d",sep="")]]$thrQ
mod_1=Mod[[paste(k,"_d",sep="")]]
#get the element for each step
ind_f=list(desc="")
nam=c("Gar","Man","Han","Hos","San")
for (i in 1:5) {
  ind_f[[paste("ind_",nam[i],1,sep="")]]=which(datath[,i]<=Qth[i])
  ind_f[[paste("N_",nam[i],1,sep="")]]=length(ind_f[[paste("ind_",nam[i],1,sep="")]])
  ind_f[[paste("ind_",nam[i],2,sep="")]]=which(datath[,i]>Qth[i])
  ind_f[[paste("N_",nam[i],2,sep="")]]=length(ind_f[[paste("ind_",nam[i],2,sep="")]])
  
}
exp1=c(3.5,2,2.5,2.5,1.5)

time_elapsed=0
ptm = proc.time()
numCores <- detectCores()-1
registerDoParallel(numCores)
for (i in 1:5) {
  res1=foreach (j=1:N2, .packages=c("stats")) %dopar% {
    beta=mod_1[[paste("beta_",nam[i],sep = "")]]
    phi=mod_1[[paste("phi_",nam[i],sep = "")]]
    sim=rep(NA,dim(beta)[1])
    ind_t=which(ind_f[[paste("ind_",nam[i],1,sep="")]]==j)
    if(length(ind_t)>0){
      mu=beta[,1]+beta[,2]*datath[j,i]^exp1[i]+beta[,3]*covP[j,i]
      sigma = phi[,1]+phi[,2]*datath[j,i]+phi[,3]*covP[j,i]
    }else{
      mu=beta[,4]+beta[,5]*datath[j,i]^exp1[i]+beta[,6]*covP[j,i]
      sigma = phi[,4]+phi[,5]*datath[j,i]+phi[,6]*covP[j,i]
    }
    
    shape = mu^2/sigma^2
    rate = mu/sigma^2
    for (k in 1:N) {
      sim[k]=rgamma(1,shape = shape[k],rate = rate[k])
    }
    sim
  }
  Q_F[[paste("Q_",k,"d",sep="")]][,,i]=t(matrix(unlist(res1), ncol = N, byrow = T))
  
}
stopImplicitCluster()


##4 day lead

k=4
datath=Q_F[[paste("VIC_For_",k,"d",sep="")]][,4:8]
tmp1=Q_F[[paste("VIC_Sim_",k,"d",sep="")]][,5:9]
covP=tmp1
Qth=Mod[[paste(k,"_d",sep="")]]$thrQ
mod_1=Mod[[paste(k,"_d",sep="")]]
#get the element for each step
ind_f=list(desc="")
nam=c("Gar","Man","Han","Hos","San")
for (i in 1:5) {
  ind_f[[paste("ind_",nam[i],1,sep="")]]=which(datath[,i]<=Qth[i])
  ind_f[[paste("N_",nam[i],1,sep="")]]=length(ind_f[[paste("ind_",nam[i],1,sep="")]])
  ind_f[[paste("ind_",nam[i],2,sep="")]]=which(datath[,i]>Qth[i])
  ind_f[[paste("N_",nam[i],2,sep="")]]=length(ind_f[[paste("ind_",nam[i],2,sep="")]])
  
}
exp1=c(2,1.5,1.5,1.5,1.5)

time_elapsed=0
ptm = proc.time()
numCores <- detectCores()-1
registerDoParallel(numCores)
for (i in 1:5) {
  res1=foreach (j=1:N2, .packages=c("stats")) %dopar% {
    beta=mod_1[[paste("beta_",nam[i],sep = "")]]
    phi=mod_1[[paste("phi_",nam[i],sep = "")]]
    sim=rep(NA,dim(beta)[1])
    ind_t=which(ind_f[[paste("ind_",nam[i],1,sep="")]]==j)
    if(length(ind_t)>0){
      mu=beta[,1]+beta[,2]*datath[j,i]^exp1[i]+beta[,3]*covP[j,i]
      sigma = phi[,1]+phi[,2]*datath[j,i]+phi[,3]*covP[j,i]
    }else{
      mu=beta[,4]+beta[,5]*datath[j,i]^exp1[i]+beta[,6]*covP[j,i]
      sigma = phi[,4]+phi[,5]*datath[j,i]+phi[,6]*covP[j,i]
    }
    
    shape = mu^2/sigma^2
    rate = mu/sigma^2
    for (k in 1:N) {
      sim[k]=rgamma(1,shape = shape[k],rate = rate[k])
    }
    sim
  }
  Q_F[[paste("Q_",k,"d",sep="")]][,,i]=t(matrix(unlist(res1), ncol = N, byrow = T))
  
}
stopImplicitCluster()


##5 day lead

k=5
datath=Q_F[[paste("VIC_For_",k,"d",sep="")]][,4:8]
tmp1=Q_F[[paste("VIC_Sim_",k,"d",sep="")]][,5:9]
covP=tmp1
Qth=Mod[[paste(k,"_d",sep="")]]$thrQ
mod_1=Mod[[paste(k,"_d",sep="")]]
#get the element for each step
ind_f=list(desc="")
nam=c("Gar","Man","Han","Hos","San")
for (i in 1:5) {
  ind_f[[paste("ind_",nam[i],1,sep="")]]=which(datath[,i]<=Qth[i])
  ind_f[[paste("N_",nam[i],1,sep="")]]=length(ind_f[[paste("ind_",nam[i],1,sep="")]])
  ind_f[[paste("ind_",nam[i],2,sep="")]]=which(datath[,i]>Qth[i])
  ind_f[[paste("N_",nam[i],2,sep="")]]=length(ind_f[[paste("ind_",nam[i],2,sep="")]])
  
}
exp1=c(2,1.5,1.5,1.5,3.5)

time_elapsed=0
ptm = proc.time()
numCores <- detectCores()-1
registerDoParallel(numCores)
for (i in 1:5) {
  res1=foreach (j=1:N2, .packages=c("stats")) %dopar% {
    beta=mod_1[[paste("beta_",nam[i],sep = "")]]
    phi=mod_1[[paste("phi_",nam[i],sep = "")]]
    sim=rep(NA,dim(beta)[1])
    ind_t=which(ind_f[[paste("ind_",nam[i],1,sep="")]]==j)
    if(length(ind_t)>0){
      mu=beta[,1]+beta[,2]*datath[j,i]^exp1[i]+beta[,3]*covP[j,i]
      sigma = phi[,1]+phi[,2]*datath[j,i]+phi[,3]*covP[j,i]
    }else{
      mu=beta[,4]+beta[,5]*datath[j,i]^exp1[i]+beta[,6]*covP[j,i]
      sigma = phi[,4]+phi[,5]*datath[j,i]+phi[,6]*covP[j,i]
    }
    
    shape = mu^2/sigma^2
    rate = mu/sigma^2
    for (k in 1:N) {
      sim[k]=rgamma(1,shape = shape[k],rate = rate[k])
    }
    sim
  }
  Q_F[[paste("Q_",k,"d",sep="")]][,,i]=t(matrix(unlist(res1), ncol = N, byrow = T))
  
}
stopImplicitCluster()

# rating curve Handia -----------------------------------------------------

i=3
exp=8/3
Q=data$data_2016$Handia_stage_10E3cms
H=data$data_2016$Handia_stage_m-data$L_Q_zero
H1=H^exp
modQ=lm(Q~H1-1)
x=seq(range(H)[1],range(H)[2],0.1)
H_sim=modQ$coefficients*x^exp
par(mfrow=c(1,1))
par(mar = c(2.5, 3, 1.5, 0.2))
temp=modQ$coefficients
limp=range(Q)
plot(data$data_2016$Handia_stage_m,Q,axes=FALSE,ann=FALSE,col=NA)
box()
grid()
axis(2,mgp=c(0, 0.5, 0))
axis(1,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow, Q"~"("*10^3*m^3*s^-1*")"), side = 2, line = 1.5, cex = 1.1,font = 1)
mtext("Water level, H (m)", side = 1, line = 1.4, cex = 1.1,font = 1)
text(x=260, y=limp[2]*0.94, labels=bquote(~ Q == .(round(temp,2))*"(H-"*H[0]*")"^.(round(exp,2))), pos=4, col="gray10",cex=1.1)
text(x=260, y=limp[2]*0.84, labels=bquote(~ H[0] == .(data$L_Q_zero)*"(m)"), pos=4, col="gray10",cex=1.1)

points(data$data_2016$Handia_stage_m,Q,pch = 21,bg="gray90",cex=1.5)
lines(x+data$L_Q_zero,H_sim,lwd=2,lty=2)
title(paste("Rating Curve at ",colnames(Q_F$Qthreshold)[i],"for 2016"))
p1.base <- recordPlot()
invisible(dev.off())
setwd(dir_JPG)
jpeg(paste("rating_curve_",colnames(Q_F$Qthreshold)[i],".jpg",sep=""), width = 6, height = 4, units = 'in', res = 300)
p1.base
dev.off()



# Plots for Garudeshwar ---------------------------------------------------

cols=c("#FD8C3B","#F16913","#D84701","#A53603","#7F2603")
cols2=c("#B8D4EE","#81B2DF","#317ABD","#2969A3","#1F4E79")
j=1
k=1
data1=Q_F[[paste("VIC_For_",k,"d",sep="")]]
ind1=which(data1$Year==Year)
par(mfrow = c(3, 1))
par(mar = c(5.8, 3.7, 0.3, 0.3))
datet=data1[ind1,1:3]
indx=which(datet$MONTH==7 & datet$DAY==1)
n2=c()
for (i in 1:length(indx)) {
  if (i==1){
    n2=c(n2,seq(indx[i],(indx[i]+61),5))
  }else{
    n2=c(n2,seq(indx[i],(indx[i]+61),5)[2:5])
  }
  
  
}

n1=c()
for (i in 1:length(ind1)) {
  tmp1=paste("0",datet$MONTH[i],sep="")
  tmp2=paste("0",datet$DAY[i],sep="")
  if (datet$MONTH[i]<10 & datet$DAY[i]<10){
    
    n1=c(n1,paste(datet$Year[i],tmp1,tmp2,sep="-"))
  }else if (datet$MONTH[i]<10) {
    n1=c(n1,paste(datet$Year[i],tmp1,datet$DAY[i],sep="-"))
  }else if (datet$DAY[i]<10){
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],tmp2,sep="-"))
  }else{
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],datet$DAY[i],sep="-"))
  }
  
  
}
limp2=c(0,0)

for (k in 1:3) {
  tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
  limp2=range(limp2,tmp_CI)
}

k=1
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(a) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()

k=2
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(b) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()

k=3
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(c) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()





p1.base <- recordPlot()
invisible(dev.off())
setwd(dir_JPG)
jpeg(paste("Forecast_",colnames(data1)[j+3],"_1_3lt.jpg",sep=""), width = 6.5, height = 7.5, units = 'in', res = 300)
p1.base
dev.off()



# Plots for Mandleshwar ---------------------------------------------------

j=2
k=1
data1=Q_F[[paste("VIC_For_",k,"d",sep="")]]
ind1=which(data1$Year==Year)
par(mfrow = c(3, 1))
par(mar = c(5.8, 3.7, 0.3, 0.3))
datet=data1[ind1,1:3]
indx=which(datet$MONTH==7 & datet$DAY==1)
n2=c()
for (i in 1:length(indx)) {
  if (i==1){
    n2=c(n2,seq(indx[i],(indx[i]+61),5))
  }else{
    n2=c(n2,seq(indx[i],(indx[i]+61),5)[2:5])
  }
  
  
}

n1=c()
for (i in 1:length(ind1)) {
  tmp1=paste("0",datet$MONTH[i],sep="")
  tmp2=paste("0",datet$DAY[i],sep="")
  if (datet$MONTH[i]<10 & datet$DAY[i]<10){
    
    n1=c(n1,paste(datet$Year[i],tmp1,tmp2,sep="-"))
  }else if (datet$MONTH[i]<10) {
    n1=c(n1,paste(datet$Year[i],tmp1,datet$DAY[i],sep="-"))
  }else if (datet$DAY[i]<10){
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],tmp2,sep="-"))
  }else{
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],datet$DAY[i],sep="-"))
  }
  
  
}
limp2=c(0,0)

for (k in 1:3) {
  tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
  limp2=range(limp2,tmp_CI)
}

k=1
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(a) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()

k=2
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(b) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()

k=3
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(c) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()





p1.base <- recordPlot()
invisible(dev.off())
setwd(dir_JPG)
jpeg(paste("Forecast_",colnames(data1)[j+3],"_1_3lt.jpg",sep=""), width = 6.5, height = 7.5, units = 'in', res = 300)
p1.base
dev.off()


# Plots for Handia --------------------------------------------------------

H=data$data_2021$Handia_stage_m-data$L_Q_zero
Q_h=modQ$coefficients*H^exp

j=3
k=1
data1=Q_F[[paste("VIC_For_",k,"d",sep="")]]
ind1=which(data1$Year==Year)
par(mfrow = c(3, 1))
par(mar = c(5.8, 3.7, 0.3, 0.3))
datet=data1[ind1,1:3]
indx=which(datet$MONTH==7 & datet$DAY==1)
n2=c()
for (i in 1:length(indx)) {
  if (i==1){
    n2=c(n2,seq(indx[i],(indx[i]+61),5))
  }else{
    n2=c(n2,seq(indx[i],(indx[i]+61),5)[2:5])
  }
  
  
}

n1=c()
for (i in 1:length(ind1)) {
  tmp1=paste("0",datet$MONTH[i],sep="")
  tmp2=paste("0",datet$DAY[i],sep="")
  if (datet$MONTH[i]<10 & datet$DAY[i]<10){
    
    n1=c(n1,paste(datet$Year[i],tmp1,tmp2,sep="-"))
  }else if (datet$MONTH[i]<10) {
    n1=c(n1,paste(datet$Year[i],tmp1,datet$DAY[i],sep="-"))
  }else if (datet$DAY[i]<10){
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],tmp2,sep="-"))
  }else{
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],datet$DAY[i],sep="-"))
  }
  
  
}
limp2=c(0,0)

for (k in 1:3) {
  tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
  limp2=range(limp2,tmp_CI)
}

k=1
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,Q_h,lty=1,lwd=2,col="gray10")
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(a) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
bias_t=round((mean(tmp)-mean(Q_h))/mean(Q_h),2)
temp=round(cor.test(tmp,Q_h,na.action("na.omit"),method="pearson")$estimate,2)
text(x=55, y=limp2[2]*0.976, labels=bquote(~ R == .(temp)), pos=4, col="gray10",cex=1.1)
text(x=53.2, y=limp2[2]*0.88, labels=bquote(~ BIAS == .(bias_t)*"%"), pos=4, col="gray10",cex=1.1)
box()

k=2
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,Q_h,lty=1,lwd=2,col="gray10")
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(b) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
bias_t=round((mean(tmp)-mean(Q_h))/mean(Q_h),2)
temp=round(cor.test(tmp,Q_h,na.action("na.omit"),method="pearson")$estimate,2)
text(x=55, y=limp2[2]*0.976, labels=bquote(~ R == .(temp)), pos=4, col="gray10",cex=1.1)
text(x=53.2, y=limp2[2]*0.88, labels=bquote(~ BIAS == .(bias_t)*"%"), pos=4, col="gray10",cex=1.1)
box()

k=3
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,Q_h,lty=1,lwd=2,col="gray10")
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(c) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
bias_t=round((mean(tmp)-mean(Q_h))/mean(Q_h),2)
temp=round(cor.test(tmp,Q_h,na.action("na.omit"),method="pearson")$estimate,2)
text(x=55, y=limp2[2]*0.976, labels=bquote(~ R == .(temp)), pos=4, col="gray10",cex=1.1)
text(x=53.2, y=limp2[2]*0.88, labels=bquote(~ BIAS == .(bias_t)*"%"), pos=4, col="gray10",cex=1.1)
box()





p1.base <- recordPlot()
invisible(dev.off())
setwd(dir_JPG)
jpeg(paste("Forecast_",colnames(data1)[j+3],"_1_3lt.jpg",sep=""), width = 6.5, height = 7.5, units = 'in', res = 300)
p1.base
dev.off()

p1.base


# Plots for Hoshangabad ---------------------------------------------------

j=4
k=1
data1=Q_F[[paste("VIC_For_",k,"d",sep="")]]
ind1=which(data1$Year==Year)
par(mfrow = c(3, 1))
par(mar = c(5.8, 3.7, 0.3, 0.3))
datet=data1[ind1,1:3]
indx=which(datet$MONTH==7 & datet$DAY==1)
n2=c()
for (i in 1:length(indx)) {
  if (i==1){
    n2=c(n2,seq(indx[i],(indx[i]+61),5))
  }else{
    n2=c(n2,seq(indx[i],(indx[i]+61),5)[2:5])
  }
  
  
}

n1=c()
for (i in 1:length(ind1)) {
  tmp1=paste("0",datet$MONTH[i],sep="")
  tmp2=paste("0",datet$DAY[i],sep="")
  if (datet$MONTH[i]<10 & datet$DAY[i]<10){
    
    n1=c(n1,paste(datet$Year[i],tmp1,tmp2,sep="-"))
  }else if (datet$MONTH[i]<10) {
    n1=c(n1,paste(datet$Year[i],tmp1,datet$DAY[i],sep="-"))
  }else if (datet$DAY[i]<10){
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],tmp2,sep="-"))
  }else{
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],datet$DAY[i],sep="-"))
  }
  
  
}
limp2=c(0,0)

for (k in 1:3) {
  tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
  limp2=range(limp2,tmp_CI)
}

k=1
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(a) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()

k=2
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(b) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()

k=3
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(c) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()





p1.base <- recordPlot()
invisible(dev.off())
setwd(dir_JPG)
jpeg(paste("Forecast_",colnames(data1)[j+3],"_1_3lt.jpg",sep=""), width = 6.5, height = 7.5, units = 'in', res = 300)
p1.base
dev.off()


# PLots for Sandiya -------------------------------------------------------

j=5
k=1
data1=Q_F[[paste("VIC_For_",k,"d",sep="")]]
ind1=which(data1$Year==Year)
par(mfrow = c(3, 1))
par(mar = c(5.8, 3.7, 0.3, 0.3))
datet=data1[ind1,1:3]
indx=which(datet$MONTH==7 & datet$DAY==1)
n2=c()
for (i in 1:length(indx)) {
  if (i==1){
    n2=c(n2,seq(indx[i],(indx[i]+61),5))
  }else{
    n2=c(n2,seq(indx[i],(indx[i]+61),5)[2:5])
  }
  
  
}

n1=c()
for (i in 1:length(ind1)) {
  tmp1=paste("0",datet$MONTH[i],sep="")
  tmp2=paste("0",datet$DAY[i],sep="")
  if (datet$MONTH[i]<10 & datet$DAY[i]<10){
    
    n1=c(n1,paste(datet$Year[i],tmp1,tmp2,sep="-"))
  }else if (datet$MONTH[i]<10) {
    n1=c(n1,paste(datet$Year[i],tmp1,datet$DAY[i],sep="-"))
  }else if (datet$DAY[i]<10){
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],tmp2,sep="-"))
  }else{
    n1=c(n1,paste(datet$Year[i],datet$MONTH[i],datet$DAY[i],sep="-"))
  }
  
  
}
limp2=c(0,0)

for (k in 1:3) {
  tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
  limp2=range(limp2,tmp_CI)
}

k=1
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(a) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()

k=2
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(b) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()

k=3
tmpV=Q_F[[paste("VIC_For_",k,"d",sep="")]][ind1,j+3]
tmp_CI=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,function(x){quantile(x,p=c(0.025,0.25,0.75,0.975))})
tmp=apply(Q_F[[paste("Q_",k,"d",sep="")]][,ind1,j], 2,mean)
# limp2=range(tmp_CI)
z1=1:length(ind1)
plot(z1,z1*0.05,type="l",lty=1,lwd=2,col="white",xlab="",ylab="",cex=1.3,cex.axis=1.1,cex.lab=1.1,axes=FALSE,ann=FALSE,ylim=limp2,xlim=c(3,60))
abline(h=0,lty=1)
polygon(c(z1,rev(z1)),c(tmp_CI[1,],rev(tmp_CI[4,])),col=col2alpha(cols[1],alpha=0.2),border=NA)
polygon(c(z1,rev(z1)),c(tmp_CI[2,],rev(tmp_CI[3,])),col=col2alpha(cols[1],alpha=0.5),border=NA)
lines(z1,tmpV,lty=1,lwd=2,col=cols2[2])
lines(z1,tmp,lty=1,lwd=2,col=cols[3])
axis(1,at=n2,labels=n1[n2],cex.axis=1.2,las=2, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
rug(x = seq(1,62,by=1), ticksize = -0.015, side = 1)
mtext(bquote("Streamflow"~"("*10^3*m^3*s^-1*")"), side = 2, line = 2, cex = 0.9,font = 1)
text(x=0.6, y=limp2[2]*0.976, labels=paste("(c) ",k,"-day lead time",sep=""), pos=4, col="black",cex=1.4,font = 2)
box()





p1.base <- recordPlot()
invisible(dev.off())
setwd(dir_JPG)
jpeg(paste("Forecast_",colnames(data1)[j+3],"_1_3lt.jpg",sep=""), width = 6.5, height = 7.5, units = 'in', res = 300)
p1.base
dev.off()


# Save the BHMC forecast for 2021 -----------------------------------------

setwd(mainDir)

if (length(Year)==1){
  list.save(Q_F, paste('data_BHM_forecast_',Year,'.rds',sep=""))
}else{
  list.save(Q_F, paste('data_BHM_forecast_',Year[1],'_',Year[length(Year)],'.rds',sep=""))
}

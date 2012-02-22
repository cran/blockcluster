#binary models
data(binarydata)
newstrategy<-cocluststrategy(nbxem=5,nbtry=2)
out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
out<-cocluster(binarydata,datatype="binary",model="pi_rho_epsilonkl",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
out<-cocluster(binarydata,datatype="binary",model="pi_rho_epsilon",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
out<-cocluster(binarydata,datatype="binary",model="pik_rhol_epsilon",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)

#gaussian models
data(gaussiandata)
newstrategy<-cocluststrategy(nbxem=5,nbtry=2)
out<-cocluster(gaussiandata,datatype="continuous",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
out<-cocluster(gaussiandata,datatype="continuous",model="pi_rho_sigma2kl",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
out<-cocluster(gaussiandata,datatype="continuous",model="pi_rho_sigma2",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
out<-cocluster(gaussiandata,datatype="continuous",model="pik_rhol_sigma2",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)

#contingency models
data(contingencydataunknown)
newstrategy<-cocluststrategy(nbxem=5,nbtry=2)
out<-cocluster(contingencydataunknown,datatype="contingency",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
out<-cocluster(contingencydataunknown,datatype="contingency",model="pi_rho_unknown",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
data(contingencydatalist)
out<-cocluster(contingencydatalist,datatype="contingency",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
out<-cocluster(contingencydatalist,datatype="contingency",model="pi_rho_known",nbcocluster=c(2,3),strategy=newstrategy,openmp=TRUE)
#binary models
data(binarydata)
newstrategy<-cocluststrategy(nbxem=5,nbtry=2)
out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3),strategy=newstrategy)
out<-cocluster(binarydata,datatype="binary",model="pi_rho_epsilonkl",nbcocluster=c(2,3),strategy=newstrategy)
out<-cocluster(binarydata,datatype="binary",model="pi_rho_epsilon",nbcocluster=c(2,3),strategy=newstrategy)
out<-cocluster(binarydata,datatype="binary",model="pik_rhol_epsilon",nbcocluster=c(2,3),strategy=newstrategy)

#gaussian models
data(gaussiandata)
newstrategy<-cocluststrategy(nbxem=5,nbtry=2)
out<-cocluster(gaussiandata,datatype="continuous",nbcocluster=c(2,3),strategy=newstrategy)
out<-cocluster(gaussiandata,datatype="continuous",model="pi_rho_sigma2kl",nbcocluster=c(2,3),strategy=newstrategy)
out<-cocluster(gaussiandata,datatype="continuous",model="pi_rho_sigma2",nbcocluster=c(2,3),strategy=newstrategy)
out<-cocluster(gaussiandata,datatype="continuous",model="pik_rhol_sigma2",nbcocluster=c(2,3),strategy=newstrategy)

#contingency models
data(contingencydataunknown)
newstrategy<-cocluststrategy(nbxem=5,nbtry=2)
out<-cocluster(contingencydataunknown,datatype="contingency",nbcocluster=c(2,3),strategy=newstrategy)
out<-cocluster(contingencydataunknown,datatype="contingency",model="pi_rho_unknown",nbcocluster=c(2,3),strategy=newstrategy)
data(contingencydatalist)
out<-cocluster(contingencydatalist,datatype="contingency",nbcocluster=c(2,3),strategy=newstrategy)
out<-cocluster(contingencydatalist,datatype="contingency",model="pi_rho_known",nbcocluster=c(2,3),strategy=newstrategy)
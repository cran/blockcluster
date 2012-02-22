#binary ss epskl
x=read.csv('/home/modal-laptop/matlabcodes/Programmes/ssbinarydataepskl.csv',head = FALSE)
data=as.matrix(x[2:dimention[1],2:dimention[2]])
dimention = dim(x)
rowlab=x[2:dimention[1],1]
collab=as.vector(t(x[1,2:dimention[2]]))
out<-cocluster(data,rowlabels=rowlab,collabels=collab,semisupervised=TRUE,datatype="binary",nbcocluster=c(2,3))
out<-cocluster(data,rowlabels=rowlab,semisupervised=TRUE,datatype="binary",nbcocluster=c(2,3)) #only row labels
out<-cocluster(data,collabels=collab,semisupervised=TRUE,datatype="binary",nbcocluster=c(2,3)) #only column labels

#binary ss eps
x=read.csv('/home/modal-laptop/matlabcodes/Programmes/ssbinarydataeps.csv',head = FALSE)
data=as.matrix(x[2:dimention[1],2:dimention[2]])
dimention = dim(x)
rowlabels=x[2:dimention[1],1]
collabels=as.vector(t(x[1,2:dimention[2]]))
out<-cocluster(data,rowlabels=rowlabels,collabels=collabels,semisupervised=TRUE,datatype="binary",model="pik_rhol_epsilon",
               nbcocluster=c(2,3))

#continuous ss sigma2kl
x=read.csv('/home/modal-laptop/matlabcodes/Programmes/ssgaussiandatasigma2kl.csv',head = FALSE)
data=as.matrix(x[2:dimention[1],2:dimention[2]])
dimention = dim(x)
rowlabels=x[2:dimention[1],1]
collabels=as.vector(t(x[1,2:dimention[2]]))
out<-cocluster(data,rowlabels=rowlabels,collabels=collabels,semisupervised=TRUE,datatype="continuous",nbcocluster=c(2,3))

#continuous ss sigma2
x=read.csv('/home/modal-laptop/matlabcodes/Programmes/ssgaussiandatasigma2.csv',head = FALSE)
data=as.matrix(x[2:dimention[1],2:dimention[2]])
dimention = dim(x)
rowlabels=x[2:dimention[1],1]
collabels=as.vector(t(x[1,2:dimention[2]]))
out<-cocluster(data,rowlabels=rowlabels,collabels=collabels,semisupervised=TRUE,datatype="continuous",nbcocluster=c(2,3))

#contingency ss unknown
x=read.csv('/home/modal-laptop/matlabcodes/Programmes/sspoissondataunknown.csv',head = FALSE)
data=as.matrix(x[2:dimention[1],2:dimention[2]])
dimention = dim(x)
rowlabels=x[2:dimention[1],1]
collabels=as.vector(t(x[1,2:dimention[2]]))
out<-cocluster(data,rowlabels=rowlabels,collabels=collabels,semisupervised=TRUE,datatype="contingency",nbcocluster=c(2,3))

#contingency ss known
x=read.csv('/home/modal-laptop/matlabcodes/Programmes/sspoissondataknown.csv',head = FALSE)
mui=read.csv('/home/modal-laptop/matlabcodes/Programmes/ssmui.dat',head = FALSE)
nuj=read.csv('/home/modal-laptop/matlabcodes/Programmes/ssnuj.dat',head = FALSE)
data=as.matrix(x[2:dimention[1],2:dimention[2]])
dimention = dim(x)
rowlabels=x[2:dimention[1],1]
collabels=as.vector(t(x[1,2:dimention[2]]))
data1=list(data,mui,nuj);
out<-cocluster(data,rowlabels=rowlabels,collabels=collabels,semisupervised=TRUE,datatype="contingency",nbcocluster=c(2,3))

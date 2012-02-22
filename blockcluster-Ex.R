pkgname <- "blockcluster"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('blockcluster')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("blockcluster")
### * blockcluster

flush(stderr()); flush(stdout())

### Name: blockcluster
### Title: Co-Clustering Package
### Aliases: blockcluster

### ** Examples

# Simple example with simulated binary data
#load data
data(binarydata)
#usage of cocluster function in its most simplest form
out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3))
#Summarize the output results
summary(out)
#Plot the original and co-clustered data
plot(out)



cleanEx()
nameEx("cocluster")
### * cocluster

flush(stderr()); flush(stdout())

### Name: cocluster
### Title: Co-Clustering function.
### Aliases: cocluster

### ** Examples

# Simple example with simulated binary data
#load data
data(binarydata)
#usage of cocluster function in its most simplest form
out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3))
#Summarize the output results
summary(out)
#Plot the original and co-clustered data
plot(out)



cleanEx()
nameEx("cocluststrategy")
### * cocluststrategy

flush(stderr()); flush(stdout())

### Name: cocluststrategy
### Title: Strategy function
### Aliases: cocluststrategy

### ** Examples

#Default strategy values

strategy<-cocluststrategy()
summary(strategy)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

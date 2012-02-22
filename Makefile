SRC_DIR = src/
COCLUSTSRCDIR = $(SRC_DIR)coclust/src/
COCLUSTLIBDIR = $(SRC_DIR)coclust/
ALGODIR = $(COCLUSTSRCDIR)Algorithms
IPDIR = $(COCLUSTSRCDIR)InputParameters
MODELDIR = $(COCLUSTSRCDIR)Models
STRATDIR = $(COCLUSTSRCDIR)Strategy
FACADEDIR = $(COCLUSTSRCDIR)CoClustFacade
COCLUSTLIB = $(LIB_DIR)libcoclust.a
RCOCLUSTLIB = $(SRC_DIR)blockcluster.so

all:
	(cd ../ && R CMD INSTALL blockcluster)

clean:	
	(cd $(SRC_DIR) && $(RM) *.o)
	(cd $(ALGODIR) && $(RM) *.o)
	(cd $(IPDIR) && $(RM) *.o)
	(cd $(MODELDIR) && $(RM) *.o)
	(cd $(FACADEDIR) && $(RM) *.o)
	(cd $(STRATDIR) && $(RM) *.o)
	(cd $(COCLUSTLIBDIR) && $(RM) *.a)
	$(RM) $(RCOCLUSTLIB)
	
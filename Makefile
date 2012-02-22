export MYCXXFLAGS = -DEIGENCONTAINERS -DRPACKAGE -fPIC
export R_HOME = /usr/lib/R
export R_SHARE_DIR = /usr/share/R/share

include $(R_HOME)/etc/Makeconf

SRC_DIR = src/
LIB_DIR = src/lib/
COCLUSTSRCDIR = $(SRC_DIR)coclust/src/
ALGODIR = $(COCLUSTSRCDIR)Algorithms
COCLUSTDIR = $(COCLUSTSRCDIR)CoClustlibrary
STKPPDIR = src/stkpp/
IPDIR = $(COCLUSTSRCDIR)InputParameters
MODELDIR = $(COCLUSTSRCDIR)Models
RCPP_PATH = /home/parmeet/R/x86_64-pc-linux-gnu-library/2.13/Rcpp/include/

COCLUSTLIB = $(LIB_DIR)libcoclust.a
STKPPLIB = $(STKPPDIR)lib/libSTKpp.a
RCOCLUSTLIB = $(SRC_DIR)Rcoclust.so

PKG_FLAGS = -I"${R_HOME}/include" -I"${RCPP_PATH}" -I"${SRC_DIR}" -I"${R_SHARE_DIR}/../include"
PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` ${COCLUSTLIB} ${STKPPLIB}

algo:
	(cd $(ALGODIR) && $(MAKE) -f algo.mk)
	
coclust:
	(cd $(COCLUSTDIR) && $(MAKE) -f coclust.mk)
	
ip:
	(cd $(IPDIR) && $(MAKE) -f ip.mk)
	
model:
	(cd $(MODELDIR) && $(MAKE) -f model.mk)
	
stkpp:
	(cd $(STKPPDIR) && $(MAKE) all)
	
CPP_SRCS = $(wildcard *.cpp)
	
CPP_OBJS = $(CPP_SRCS:%.cpp=$(SRC_DIR)%.o)

$(SRC_DIR)%.o: $(SRC_DIR)%.cpp
	$(CXX) $(PKG_FLAGS) $(MYCXXFLAGS) -Wall -pedantic $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) $< -c -o $@
	
.PHONY: algo coclust ip model stkpp all 

all: algo coclust ip model stkpp $(RCOCLUSTLIB)

$(RCOCLUSTLIB): $(CPP_OBJS) $(COCLUSTLIB)
	$(CXX) $(SHLIB_CXXLDFLAGS) -o $@ $(CPP_OBJS) $(PKG_LIBS)

clean:	
	(cd $(SRC_DIR) && $(RM) *.o)
	(cd $(ALGODIR) && $(RM) *.o)
	(cd $(IPDIR) && $(RM) *.o)
	(cd $(MODELDIR) && $(RM) *.o)
	(cd $(COCLUSTDIR) && $(RM) *.o)
	(cd $(COCLUSTSRCDIR) && $(RM) *.a)
	(cd $(STKPPDIR) && $(MAKE) clean)
	$(RM) $(COCLUSTLIB) $(RCOCLUSTLIB)
	
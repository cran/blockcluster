include $(R_HOME)/etc${R_ARCH}/Makeconf
MYCXXFLAGS = -DNDEBUG -DRPACKAGE
LIB = libcoclust.a
SRCS = $(wildcard src/Algorithms/*.cpp\
                  src/CoClustFacade/*.cpp\
                  src/Strategy/*.cpp\
                  src/InputParameters/*.cpp\
                  src/Models/*.cpp)

OBJS= $(SRCS:%.cpp=%.o)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) ${CPICFLAGS} $(MYCXXFLAGS) $< -c -o $@
	
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rc $@ $?
  
.PHONY: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o

#-----------------------------------------------------------------------
#     Copyright (C) 2004-2007  Serge Iovleff
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : Serge.Iovleff@stkpp.org
#
#-----------------------------------------------------------------------
#
# Purpose:  Main makefile of the Project STKpp.
# Author:   Serge Iovleff, serge.iovleff@stkpp.org
#
#-----------------------------------------------------------------------
#
.PHONY: all \
        obj \
        lib \
        stkernel \
        dmanager \
        aamodels \
        regress \
        reduct \
        algebra \
        statistik \
        libclean \
        objclean \
        testclean \
        clean

#-----------------------------------------------------------------------
# Macros and definitions
include ./makevars

#-----------------------------------------------------------------------
# main Paths
PROJECTS_DIR = ./projects/
GUI_DIR = ./gui/

#-----------------------------------------------------------------------
# Overide Paths BIN and LIB
BIN_DIR = ./bin/
LIB_DIR = ./

#-----------------------------------------------------------------------
# specific projects paths
STKERNEL_DIR = $(PROJECTS_DIR)STKernel/
DMANAGER_DIR = $(PROJECTS_DIR)DManager/
ALGEBRA_DIR = $(PROJECTS_DIR)Algebra/
STATISTIK_DIR = $(PROJECTS_DIR)STatistiK/
AAMODELS_DIR = $(PROJECTS_DIR)AAModels/
REGRESS_DIR = $(PROJECTS_DIR)Regress/
REDUCT_DIR = $(PROJECTS_DIR)Reduct/
MODEL_DIR = $(PROJECTS_DIR)Model/
SDK_DIR = $(PROJECTS_DIR)Sdk/

SRC_STKERNEL_DIR = $(STKERNEL_DIR)src/
SRC_DMANAGER_DIR = $(DMANAGER_DIR)src/
SRC_ALGEBRA_DIR = $(ALGEBRA_DIR)src/
SRC_STATISTIK_DIR = $(STATISTIK_DIR)src/
SRC_AAMODELS_DIR = $(AAMODELS_DIR)src/
SRC_REGRESS_DIR = $(REGRESS_DIR)src/
SRC_REDUCT_DIR = $(REDUCT_DIR)src/
SRC_MODEL_DIR = $(MODEL_DIR)src/
SRC_SDK_DIR = $(SDK_DIR)src/
#-----------------------------------------------------------------------
# projects binaries
STKERNEL_OBJS = $(patsubst $(SRC_STKERNEL_DIR)%.cpp, $(BIN_DIR)%.o, $(wildcard $(SRC_STKERNEL_DIR)/*.cpp))
ALGEBRA_OBJS = $(patsubst $(SRC_ALGEBRA_DIR)%.cpp, $(BIN_DIR)%.o, $(wildcard $(SRC_ALGEBRA_DIR)/*.cpp))
DMANAGER_OBJS = $(patsubst $(SRC_DMANAGER_DIR)%.cpp, $(BIN_DIR)%.o, $(wildcard $(SRC_DMANAGER_DIR)/*.cpp))
STATISTIK_OBJS = $(patsubst $(SRC_STATISTIK_DIR)%.cpp, $(BIN_DIR)%.o, $(wildcard $(SRC_STATISTIK_DIR)/*.cpp))
REGRESS_OBJS =$(patsubst $(SRC_REGRESS_DIR)%.cpp, $(BIN_DIR)%.o, $(wildcard $(SRC_REGRESS_DIR)/*.cpp))
AAMODELS_OBJS = $(patsubst $(SRC_AAMODELS_DIR)%.cpp, $(BIN_DIR)%.o, $(wildcard $(SRC_AAMODELS_DIR)/*.cpp))
REDUCT_OBJS = $(patsubst $(SRC_REDUCT_DIR)%.cpp, $(BIN_DIR)%.o, $(wildcard $(SRC_REDUCT_DIR)/*.cpp))
MODEL_OBJS = $(patsubst $(SRC_MODEL_DIR)%.cpp, $(BIN_DIR)%.o, $(wildcard $(SRC_MODEL_DIR)/*.cpp))
SDK_OBJS = $(patsubst $(SRC_SDK_DIR)%.cpp, $(BIN_DIR)%.o, $(wildcard $(SRC_SDK_DIR)/*.cpp))
# all binaries
CPP_OBJS = \
  $(STKERNEL_OBJS) \
  $(AAMODELS_OBJS) \
  $(MODEL_OBJS) \
  $(SDK_OBJS)
  

# library name
CPP_LIB = $(LIB_DIR)$(LIB_STKPP_A)

#-----------------------------------------------------------------------
# create the library
$(CPP_LIB): $(CPP_OBJS)
	$(AR) rv $@ $?

#-----------------------------------------------------------------------
# obj compiles the objects in all projects
obj:
	cd $(SRC_STKERNEL_DIR); $(MAKE) all
	#cd $(SRC_DMANAGER_DIR); $(MAKE) all
	cd $(SRC_AAMODELS_DIR); $(MAKE) all
	cd $(SRC_MODEL_DIR); $(MAKE) all
	cd $(SRC_SDK_DIR); $(MAKE) all
	#cd $(SRC_ALGEBRA_DIR); $(MAKE) all
	#cd $(SRC_STATISTIK_DIR); $(MAKE) all
	#cd $(SRC_REGRESS_DIR); $(MAKE) all
	#cd $(SRC_REDUCT_DIR); $(MAKE) all

#-----------------------------------------------------------------------
# create the library
lib: $(CPP_LIB)

#-----------------------------------------------------------------------
# all compiles the objects and create the library stkpp.a
all: obj lib

#-----------------------------------------------------------------------
# rule if we want to compile tests programs
tests: obj lib
	cd tests/; $(MAKE) all

#-----------------------------------------------------------------------
# rules if we want to compile only one sub-project (assuming nothing
# else is modified)
stkernel:
	cd $(SRC_STKERNEL_DIR); $(MAKE) all

dmanager:
	cd $(SRC_DMANAGER_DIR); $(MAKE) all

algebra:
	cd $(SRC_ALGEBRA_DIR); $(MAKE) all

statistik:
	cd $(SRC_STATISTIK_DIR); $(MAKE) all

aamodels:
	cd $(SRC_AAMODELS_DIR); $(MAKE) all

regress:
	cd $(SRC_REGRESS_DIR); $(MAKE) all

reduct:
	cd $(SRC_REDUCT_DIR); $(MAKE) all

#-----------------------------------------------------------------------
# rules for cleaning objects
clean: libclean objclean testclean

libclean:
	cd $(LIB_DIR); rm -f *.a *.so *.dll

objclean:
	cd $(BIN_DIR); rm -f *.o

testclean:
	cd $(BIN_DIR); rm -f test*

#-----------------------------------------------------------------------
# gui project
qt:
	cd ./gui/qt/; qmake; $(MAKE)

model:
	cd ./gui/model/src; $(MAKE) all

gui: model qt

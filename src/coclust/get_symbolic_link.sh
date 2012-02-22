#!/bin/bash
##############################################
# define function to make symbolic link of files from a directory

function make_symbolic_link
{
# loop over header files
for file in $1/*.h
do 
link=./src/$(basename $1)/$(basename $file)
ln -s ../../$file $link
#echo $file
done

# loop over header files
for file in $1/*.cpp
do 
link=./src/$(basename $1)/$(basename $file)
ln -s ../../$file $link
#echo $link
done
}
make_symbolic_link ../../../coclust/coclust/src/Algorithms
make_symbolic_link ../../../coclust/coclust/src/CoClustFacade
make_symbolic_link ../../../coclust/coclust/src/enumerations
make_symbolic_link ../../../coclust/coclust/src/Initialization
make_symbolic_link ../../../coclust/coclust/src/InputParameters
make_symbolic_link ../../../coclust/coclust/src/Models
make_symbolic_link ../../../coclust/coclust/src/typedefs
make_symbolic_link ../../../coclust/coclust/src/Strategy

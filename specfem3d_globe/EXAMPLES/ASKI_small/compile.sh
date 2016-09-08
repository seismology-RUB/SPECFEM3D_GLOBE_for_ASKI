#!/bin/bash
#
# ASKI small example
#
# script compiles only!
##################################################

echo "compiling example: `date`"
currentdir=`pwd`

echo "directory: $currentdir"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

# compiles executables in root directory
# using default configuration
cd ../../

# compiles for a forward simulation
cp $currentdir/DATA/Par_file DATA/Par_file
make clean
make xmeshfem3D xspecfem3D

# backup of constants setup
cp setup/* $currentdir/OUTPUT_FILES/
cp OUTPUT_FILES/values_from_mesher.h $currentdir/OUTPUT_FILES/values_from_mesher.h.compilation
cp DATA/Par_file $currentdir/OUTPUT_FILES/

cd $currentdir

# copying the executables ist not necessary, since ./bin contains symlinks to ../../bin/xmeshfem3D and ../../bin/xspecfem3D

echo `date`


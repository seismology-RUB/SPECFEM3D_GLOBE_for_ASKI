#!/bin/bash

###########################################################
# USER PARAMETERS

## 4 CPUs
CPUs=4

###########################################################


BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# script to run the mesher and the solver
# read DATA/Par_file to get information about the run
# compute total number of nodes needed
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

if [ ! "$numnodes" == "$CPUs" ]; then
  echo "error: Par_file for $numnodes CPUs"
  exit 1
fi

rm -rf OUTPUT_FILES/*

# backup files used for this simulation
cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp ../../setup/constants.h OUTPUT_FILES/


# restoring files from mesher that were backuped at mesh generation
cp $BASEMPIDIR/*.txt OUTPUT_FILES/
cp $BASEMPIDIR/values_from_mesher.h OUTPUT_FILES/
cp $BASEMPIDIR/output_mesher.txt OUTPUT_FILES/


##
## forward simulation
##

sleep 2

echo
echo `date`
echo starting run in current directory $PWD
echo

mpirun -np $numnodes $PWD/bin/xspecfem3D

echo "finished successfully"
echo `date`


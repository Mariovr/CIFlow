#!/bin/sh
#
#
#PBS -N ciflowtest
#PBS -o output.file
#PBS -e error.file
#PBS -l walltime=00:59:59
#PBS -q debug
#PBS -m be
#PBS -l nodes=1:ppn=8
#

SRC_HOME_DIR=$VSC_HOME/devel/CIFlow/rundir
echo $PWD

cd $SRC_HOME_DIR
rm -rf $VSC_SCRATCH_NODE/ciflow
mkdir $VSC_SCRATCH_NODE/ciflow

cp ./flow.dat $VSC_SCRATCH_NODE/ciflow/flow.dat
cp ./ciflow.x $VSC_SCRATCH_NODE/ciflow/ciflow.x
cp ./output.dat $VSC_SCRATCH_NODE/ciflow/output.dat

cd $VSC_SCRATCH_NODE/ciflow
./ciflow.x < flow.dat
cd $SRC_HOME_DIR

rm -rf debug_run
mkdir debug_run
cp -r $VSC_SCRATCH_NODE/ciflow   ./debug_run/.

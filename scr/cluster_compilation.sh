#!/bin/bash
#PBS -N CIFlow compilation
#PBS -q debug

module load scripts
module load GCC/4.9.1
module load icc/2015.0.090
module load ifort/2015.0.090
module load iccifort/2015.0.090
module load impi/5.0.1.035-iccifort-2015.0.090
module load iimpi/7.1.2
module load imkl/11.2.0.090-iimpi-7.1.2
module load ictce/7.1.2
module load zlib/1.2.8-ictce-7.1.2
module load Szip/2.1-ictce-7.1.2
module load HDF5/1.8.14-ictce-7.1.2-serial
module load bzip2/1.0.6-ictce-7.1.2
module load ncurses/5.9-ictce-7.1.2
module load libreadline/6.3-ictce-7.1.2
module load Boost/1.57.0-intel-2015a-Python-2.7.9
module load Python/2.7.8-ictce-7.1.2
module load arpack-ng/3.1.5-ictce-7.1.2-mt

SRC_HOME_DIR=$VSC_HOME/devel/CIFlow
echo $PWD
BINARY_NAME="ciflow_"
CLUSTER=`qstat -q | grep 'server' | sed -e 's/server: *[a-z]*[0-9]*\.\([a-z]*\)\.gent.vsc/\1/g'`
BINARY_NAME+=$CLUSTER
BINARY_NAME+=".x"

echo "Building for cluster: $CLUSTER" 
cd $SRC_HOME_DIR
rm -rf $VSC_SCRATCH_NODE/ciflow
mkdir $VSC_SCRATCH_NODE/ciflow
cp -r ./src $VSC_SCRATCH_NODE/ciflow/.
cp -r ./include $VSC_SCRATCH_NODE/ciflow/.
cp -r ./tests $VSC_SCRATCH_NODE/ciflow/.
cp Makefile $VSC_SCRATCH_NODE/ciflow/.
cd $VSC_SCRATCH_NODE/ciflow && make program && make tests
cd $SRC_HOME_DIR && cp $VSC_SCRATCH_NODE/ciflow/bin/ciflow.x ./bin/$BINARY_NAME

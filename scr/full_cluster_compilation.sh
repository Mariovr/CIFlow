#!/bin/bash

CLUSTERS=( "raichu" "delcatty" )
for CLUSTER in ${CLUSTERS[@]} 
do
    echo $CLUSTER
    module swap cluster/$CLUSTER
    qsub cluster_compilation.sh
done
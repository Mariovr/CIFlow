#!/bin/bash
make
cd bin
shopt -s nullglob
for exe in *.x
do
    ./$exe
done
exit 0

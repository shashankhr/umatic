#!/bin/sh
#Shell script to create an input routine from a structure definition
#Some manual editing will be required subsequently
if [[ $# != 1 ]]
then
   echo $#
   echo -e "$0: Create input routine from a structure definition.\nUsage $0 filename.h "
   echo " where filename.h contains a structure definition suitable for control inputs"
   echo "Requires: compiled make_structure_reader and make_structure_output programs "
   exit
fi
n=$1
echo $n
o=`cat $n | make_structure_output | tail -1 `
i=`cat $n | make_structure_reader | tail -1 `
echo $o
echo $i

sed -e'/MANUALLY ADD the header/a\
#include "'$n'"' $o > $o.fix
mv $o.fix $o

sed -e'/MANUALLY ADD the header/a\
#include "'$n'"' $i > $i.fix
mv $i.fix $i


#!/bin/sh
#
# Shell script to take a table of values
# with the ca_matprop commands in headings
# and create a series of subdirectories
# with the substituted ca_matprop.in files
# Subdirectory names must be in the first column
#
echo " "
echo " "
echo "****************************"
echo "*      makeruns.sh         *"
echo "****************************"
echo " "
echo " "

# get the file name to read from
in=$1
# print usage if no argument
if [[ ! -n $1 ]]
then
   echo "Shell script to take a table of values"
   echo "with the ca_matprop commands in headings"
   echo "and create a series of subdirectories"
   echo "with the substituted ca_matprop.in files"
   echo "Subdirectory names must be in the first column"
   echo " "
   echo "Usage: makeruns.sh infile"
   echo " "
   echo "from a directory containing template ctrl files."
   echo "The infile should contain a table of parameters."
   echo "The first line contains:"
   echo "experimentname caparameter1 ca_parameter2 etc"
   echo "The remaining lines conatain:"
   echo "runname value1 value2 etc"
   echo "where experimentname and runname are your names and"
   echo "ca_parameters are valid option fields for ca control files"
   exit
fi

echo "In-file: " $in
declare a names
declare a vals
# exit if filename is not ordinary.
if [[ ! -f $in ]]
then
  echo "The file $in does not exist or something is wrong."
  echo "makeruns.sh with no argument for the usage message."
  echo "Exiting ... "
  exit
fi

#fix it since it was probably saved from Excel in windoze
to_unix $in $in
#flag for the header line
f=0
while read a 
do
   if [[ $f -eq 0 ]] 
then
      #it is the first line
      f=1
      echo "First line: \n" $a
      #set the array of commands; subdir names are in first column
      #set -A names $a
      names=($a)
   else
      #it is not the first line. Make a temporary copy of the template file
      cp ca_matprop.in ca_m_runs.in
      cp ca_geoplus.in ca_g_runs.in
      cp ca_ctrl.in ca_c_runs.in
      #set the value array from the read line
      vals=($a)
      echo $a
      if [[ ! -d ${vals[0]} ]] 
      then
         mkdir ${vals[0]}
      fi
      #
      #substitute ca_matprop file
      #
      #flag for name element
      j=0
      for i in ${names[*]}
      do
         if [[ $j -ne 0 ]] 
      then
            n=${names[$j]}" "
            v=${vals[$j]}" "
            #substitute this command
            sed -e"s/^$n *[0-9\.Ee+-]*/$n $v /" ca_m_runs.in > ca_m_runs.tmp
            #wrap the file to be ready for next substitution
            mv ca_m_runs.tmp ca_m_runs.in
         fi
         ((j=$j + 1))
      done
      #put it in the directory
      mv ca_m_runs.in ${vals[0]}/ca_matprop.in
      #
      #substitute ca_geoplus file
      #
      #flag for name element
      j=0
      for i in ${names[*]}
      do
         if [[ $j -ne 0 ]] 
      then
            n=${names[$j]}" "
            v=${vals[$j]}" "
            #substitute this command
            sed -e"s/^$n *[0-9\.Ee+]*/$n $v /" ca_g_runs.in > ca_g_runs.tmp
            #wrap the file to be ready for next substitution
            mv ca_g_runs.tmp ca_g_runs.in
         fi
         ((j=$j + 1))
      done
      #put it in the directory
      mv ca_g_runs.in ${vals[0]}/ca_geoplus.in
      #
      # substitute ca_ctrl file
      #
      #flag for name element
      j=0
      for i in ${names[*]}
      do
         if [[ $j -ne 0 ]] 
      then
            n=${names[$j]}" "
            v=${vals[$j]}" "
            #substitute this command
            sed -e"s/^$n *[-0-9\.\Ee+]*/$n $v /" ca_c_runs.in > ca_c_runs.tmp
            #wrap the file to be ready for next substitution
            mv ca_c_runs.tmp ca_c_runs.in
         fi
         ((j=$j + 1))
      done
      #now substitute base file name in ca_ctrl.in
      n="BaseFileName"
      v="${names[0]}_${vals[0]} \/\*base filename\*\/"
      echo $v
      sed -e"s/$n.*$/$n $v /" ca_c_runs.in > ${vals[0]}/ca_ctrl.in
      cp ca_run ${vals[0]}/ca_run
      #copy any input files for mould
      cp *.inp ${vals[0]}
      #copy any mat prop files
      cp props*.in ${vals[0]}
      rm ca_c_runs.in
   fi
#this is where the input file is given to the read command !?!?!?!
done < $in

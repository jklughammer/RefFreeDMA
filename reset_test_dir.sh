#!/bin/bash
set -e
dir=$1

shopt -s extglob

read -p "ATTENTION:THIS WILL WIPE THE INDICATED DIRECTORY. Do you really want to reset the directory $dir (y/n)? " RESP
if [ $RESP = "y" ]; then
  cd $dir 
  rm -r !(meta|unmapped_bam) 2>/dev/null
  rm *.done 2>/dev/null
  echo "$dir is now reset."
elif [ $RESP = "n" ]; then
  echo "Abort!"
else
  echo "Unknown entry.(y/n) needed.Exiting!"
fi 


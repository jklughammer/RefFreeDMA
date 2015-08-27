#!/bin/bash

dir=$1

shopt -s extglob
cd $dir
rm -r  !(meta|unmapped_bam) 2>/dev/null
rm *.done 2>/dev/null

echo "$dir is now reset."

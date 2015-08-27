#!/bin/bash

dir=RefFreeDMA_test

shopt -s extglob
cd $dir
rm -r  !(meta|unmapped_bam) 
rm *.done

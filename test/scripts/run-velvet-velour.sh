#!/bin/bash

# run both velvet and velour and compare output

#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

PREFIX=$1
KMER=$2
INPUT=$3

VELOUR=../../../velour

VELVETH=../../../../velvet/velveth
VELVETG=../../../../velvet/velvetg

COMPARE=../../scripts/compare-pregraph-simple.py

# velvet

if [ -d velvet.$1 ] ; then
    rm -rf velvet.$1
fi
mkdir velvet.$1
$VELVETH velvet.$1 $2 $3 > velvet.$1/thelog
$VELVETG velvet.$1 >> velvet.$1/thelog

# velour

if [ -d velour.$1 ] ; then
    rm -rf velour.$1
fi
mkdir velour.$1
$VELOUR velour.$1 $2 $3 > velour.$1/thelog

# compare

$COMPARE velvet.$1/PreGraph velour.$1/PreGraph


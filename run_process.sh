#!/bin/bash

sampleFile=$1
nCPU=3

cat $sampleFile | sed "s/,/ /g" | xargs -n5 -P $nCPU ~/tools/sctools/process.R

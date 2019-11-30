#!/bin/bash

#list of L parameters
L="5,7,11,13,17,19,20,40,50,80,100,120,200"

#list of distance parameters
declare -a dist=("mean" "max" "min")


#perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/1_2 2 'mean' "$L" ~/Documents/results/k2/mean/
#arguments: [location of vcf] [k] [dist] [list of L] [output folder]

for d in "${dist[@]}"
do
  perl "~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5_6 2 $d $L ~/Documents/results/k2/$d &"
  perl "~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5_6 3 $d $L ~/Documents/results/k3/$d &"
  perl "~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5_6 4 $d $L ~/Documents/results/k4/$d &"
  perl "~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5_6 5 $d $L ~/Documents/results/k5/$d &"
  wait
  echo "done $d"
done

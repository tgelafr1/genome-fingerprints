#!/bin/bash

#list of L parameters
L="5,7,11,13,17,19,20,40,50,80,100,120,200"

#list of distance parameters
# declare -a dist=("mean" "max" "min")


#perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/1_2 2 'mean' "$L" ~/Documents/results/k2/mean/
#arguments: [location of vcf] [k] [dist] [list of L] [output folder]

# mean
# perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 2 mean "$L" ~/Documents/results/k2/mean
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 3 mean "$L" ~/Documents/results/k3/mean
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 4 mean "$L" ~/Documents/results/k4/mean
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 5 mean "$L" ~/Documents/results/k5/mean
# max
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 2 max "$L" ~/Documents/results/k2/max
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 3 max "$L" ~/Documents/results/k3/max
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 4 max "$L" ~/Documents/results/k4/max
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 5 max "$L" ~/Documents/results/k5/max

# min
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 2 min "$L" ~/Documents/results/k2/min
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 3 min "$L" ~/Documents/results/k3/min
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 4 min "$L" ~/Documents/results/k4/min
perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 5 min "$L" ~/Documents/results/k5/min

# for d in "${dist[@]}"
# do
#   perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 2 "$d" "$L" ~/Documents/results/k2/"$d" &
#   perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 3 "$d" "$L" ~/Documents/results/k3/"$d" &
#   perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 4 "$d" "$L" ~/Documents/results/k4/"$d" &
#   perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/5 5 "$d" "$L" ~/Documents/results/k5/"$d" &
#   echo "finished $d "
# done

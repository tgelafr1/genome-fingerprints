#!/bin/bash

#list of folders that contain 1000 genomes multisample vcf
declare -a arr=("1_2" "3_4" "5_6" "7_8" "9_10" "11_12" "13_14" "15_16" "17_18" "19_20" "21_22")

#list of L parameters
L="5,7,11,13,17,19,20,40,50,80,100,120,200"

#perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/1_2 2 'mean' "$L" ~/Documents/results/k2/mean/
#arguments: [location of vcf] [k] [dist] [list of L] [output folder]

## dist=mean

#run for k=2, dist=mean
for i in "${arr[@]}"
do
	perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 2 'mean' "$L" ~/Documents/results/k2/mean &
done
wait
echo "\nk=2, dist=mean\n"

#run for k=3, dist=mean
for i in "${arr[@]}"
do
	perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 3 'mean' "$L" ~/Documents/results/k3/mean &
done
wait
echo "\nk=3, dist=mean\n"

#run for k=4, dist=mean
for i in "${arr[@]}"
do
        perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 4 'mean' "$L" ~/Documents/results/k4/mean &
done
wait
echo "\nk=4, dist=mean\n"

#run for k=5, dist=mean
for i in "${arr[@]}"
do
	perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 5 'mean' "$L" ~/Documents/results/k5/mean &
done
wait
echo "\nk=5, dist=mean\n"

## dist=max

#run for k=2, dist=max
for i in "${arr[@]}"
do
	perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 2 'max' "$L" ~/Documents/results/k2/max &
done
wait
echo "\nk=2, dist=max\n"

#run for k=3, dist=max
for i in "${arr[@]}"
do
        perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 3 'max' "$L" ~/Documents/results/k3/max &
done
wait
echo "\nk=3, dist=max\n"

#run for k=4, dist=max
for i in "${arr[@]}"
do
        perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 4 'max' "$L" ~/Documents/results/k4/max &
done
wait
echo "\nk=4, dist=max\n"

#run for k=5, dist=max
for i in "${arr[@]}"
do
        perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 5 'max' "$L" ~/Documents/results/k5/max &
done
wait
echo "\nk=5, dist=max\n"

## dist=min

#run for k=2, dist=min
for i in "${arr[@]}"
do
	perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 2 'min' "$L" ~/Documents/results/k2/min &
done
wait
echo "\nk=2, dist=min\n"

#run for k=3, dist=min
for i in "${arr[@]}"
do
        perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 3 'min' "$L" ~/Documents/results/k3/min &
done
wait
echo "\nk=3, dist=min\n"

#run for k=4, dist=min
for i in "${arr[@]}"
do
        perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 4 'min' "$L" ~/Documents/results/k4/min &
done
wait
echo "\nk=4, dist=min\n"

#run for k=5, dist=min
for i in "${arr[@]}"
do
        perl ~/github/genome-fingerprints/bin/computeDMF-multisample.pl ~/Documents/1000genomes/"$i" 5 'min' "$L" ~/Documents/results/k5/min &
done
wait
echo "\nk=5, dist=min\n"

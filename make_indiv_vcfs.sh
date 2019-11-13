#!/bin/bash

vcf=$1

for i in {46..60}; do
	echo $i
	cat $vcf | cut -d $'\t' -f 1-9,$i > "hg${i}.vcf"
done

#!/bin/bash
# original fingerprinting
perl bin_orig/computeDMF.pl data/results/chr1_orig ~/Documents/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
perl bin_orig/computeDMF.pl data/results/chr2_orig ~/Documents/1000genomes/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# modified fingerprinting
perl bin/computeDMF.pl data/results/chr1 ~/Documents/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
perl bin/computeDMF.pl data/results/chr2 ~/Documents/1000genomes/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

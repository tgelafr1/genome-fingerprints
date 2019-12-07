# genome-fingerprints
Software for creating and comparing genome fingerprints.
Forked from https://github.com/gglusman/genome-fingerprints
Tested on Ubuntu 18.04.3 LTS (Bionic Beaver)

## dependencies
`perl (v.5.26.2)`

`python 3`

`numpy`, `pandas`, `joblib`, `tqdm`, `matplotlib`, `seaborn`, `scikit-learn`

## 1000 Genomes data
VCF files were obtained from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

Example multisample VCF file: `ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`

## Choosing samples and populations
Header file with complete list of all 2504 individuals and their corresponding population are in `data/actual_sample_names.csv`

You can compute fingerprints on a subset of the data by using `extract_population.ipynb`

## Computing fingerprints
perl script that computes fingerprints from `*.vcf.gz` is bin/computeDMF.pl

python script that computes multiple fingerprints from a multisample VCF: `make_indiv_vcfs.py`

Example usage:

`python make_indiv_vcfs.py bin/computeDMF.pl 3 mean data/results_chr21/ data/sample_names_6_pop.csv data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 5,7,11,13,17,19,20,40,50,80,100,120,200`

## Evaluation of algorithm
Use `Analysis.ipynb` to compare fingerprints to populations, do PCA, and classify samples to populations

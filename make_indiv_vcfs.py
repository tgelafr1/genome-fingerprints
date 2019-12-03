from pathlib import Path
import os
import pandas as pd
import sys
import subprocess
import time

import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm


# read in arguments
# header = pd.read_csv(sys.argv[1])
# vcf = sys.argv[2]
# compute = sys.argv[3]

# parameters
compute = '/home/ubuntu/github/genome-fingerprints/bin/computeDMF.pl'
results_path = Path('/home/ubuntu/Documents/results_chr21')
header = pd.read_csv('/home/ubuntu/github/genome-fingerprints/data/sample_names_6_pop.csv')

vcf = '/home/ubuntu/github/genome-fingerprints/data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
k = str(5)
dist = 'mean'
format = 'vcf'
L = '5,7,11,13,17,19,20,40,50,80,100,120,200'

# create output folder
output_path = results_path / ('k' + k) / dist
if not os.path.exists(output_path):
    os.makedirs(output_path)


# get set of inputs we want our function to iterate over in parallel
names = header.Sample
inputs = enumerate(tqdm(names)) #tqdm for progress bar

# get number of cores/jobs
# num_cores = multiprocessing.cpu_count()
num_cores = 2

def individual_vcf(i, sample):

    print(sample)

    sample_vcf_file = str(output_path / sample )

    fingerprint_file = str(output_path/ sample)

    subprocess.run("zcat " + vcf + " | cut -f 1-9," + str(i+10) + " > " + sample_vcf_file, shell = True)
    subprocess.call("perl " + compute + " " + fingerprint_file + " " + sample_vcf_file + " " + k + " " + dist + " " + format + " " + L, shell = True)
    subprocess.call("rm " + sample_vcf_file, shell = True)


if __name__ == "__main__":
    # run individual_vcf in parallel
    Parallel(n_jobs=num_cores)(delayed(individual_vcf)(i[0],i[1]) for i in inputs)

"""Python script to compute fingerprints with parallelization capability

Example usage:
python make_indiv_vcfs.py bin/computeDMF.pl 3 mean results_chr21/ data/sample_names_6_pop.csv data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 5,7,11,13,17,19,20,40,50,80,100,120,200
"""

from pathlib import Path
import os
import sys
import subprocess
import time

import pandas as pd
import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm

def individual_vcf(i, sample, compute, results_path, header, vcf, k, dist, format, L):
    """Computes fingerprint on 1 individual (sample)
    """
    print(sample)

    sample_vcf_file = str(output_path / sample )

    fingerprint_file = str(output_path/ sample)

    subprocess.run("zcat " + vcf + " | cut -f 1-9," + str(i+10) + " > " + sample_vcf_file, shell = True)
    subprocess.call("perl " + compute + " " + fingerprint_file + " " + sample_vcf_file + " " + k + " " + dist + " " + 'vcf' + " " + L, shell = True)
    subprocess.call("rm " + sample_vcf_file, shell = True)


if __name__ == "__main__":

    # parameters
    # compute = 'bin/computeDMF.pl'
    # k = str(5)
    # dist = 'mean'
    # results_path = Path('results_chr21/')
    # header = pd.read_csv('data/sample_names_6_pop.csv')
    # vcf = 'data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
    # L = '5,7,11,13,17,19,20,40,50,80,100,120,200'
    # format = 'vcf'
    compute = sys.argv[1]
    k = str(sys.argv[2])
    dist = sys.argv[3]
    results_path = Path(sys.argv[4])
    header = pd.read_csv(sys.argv[5])
    vcf = sys.argv[6]
    L = sys.argv[7]

    # get number of cores/jobs
    num_cores = 2

    # create output folder
    output_path = results_path / ('k' + k) / dist
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # get set of inputs we want our function to iterate over in parallel
    names = header.Sample
    inputs = enumerate(tqdm(names)) #tqdm for progress bar

    # run individual_vcf in parallel
    Parallel(n_jobs=num_cores)(delayed(individual_vcf)(i[0],i[1],compute, results_path, header, vcf, k, dist, format, L) for i in inputs)

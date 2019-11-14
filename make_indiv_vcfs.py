import pandas as pd
import sys
import subprocess

import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm


# read in arguments
header = pd.read_csv(sys.argv[1])
vcf = sys.argv[2]
compute = sys.argv[3]

# get set of inputs we want our function to iterate over in parallel
names = header.Sample
inputs = enumerate(tqdm(names)) #tqdm for progress bar

# get number of cores/jobs
num_cores = multiprocessing.cpu_count()

def individual_vcf(i, sample):
    vcf_file = sample + ".vcf"
    subprocess.run("zcat " + vcf + " | cut -f 1-9," + str(i+10) + " > " + vcf_file, shell = True)
    subprocess.call("perl " + compute + " "  +  sample + " " + vcf_file + " sum", shell = True)
    subprocess.call("rm " + vcf_file, shell = True)
    print(sample)



if __name__ == "__main__":
    # run individual_vcf in parallel
    Parallel(n_jobs=num_cores)(delayed(individual_vcf)(i[0],i[1]) for i in inputs)

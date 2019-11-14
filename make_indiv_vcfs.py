import pandas as pd
import sys
import subprocess

header = pd.read_csv(sys.argv[1])
names = header.Sample

vcf = sys.argv[2]
compute = sys.argv[3]

for i, sample in enumerate(names):
    vcf_file = sample + ".vcf"
    subprocess.run("zcat " + vcf + " | cut -f 1-9," + str(i+10) + " > " + vcf_file, shell = True)
    subprocess.call("perl " + compute + " "  +  sample + " " + vcf_file + " sum", shell = True)
    subprocess.call("rm " + vcf_file, shell = True)

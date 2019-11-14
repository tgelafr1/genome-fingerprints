import pandas as pd
import sys
import subprocess

header = pd.read_csv(sys.argv[1])
names = header.Sample
vcf = sys.argv[2]
compute = sys.argv[3]

for i, sample in enumerate(names):
    print(sample)
    vcf_file = sample + ".vcf"
    subprocess.call(["cat", vcf, "|", "cut", "-f", "1-9", i+10, ">", vcf_file])
    subprocess.call(["perl", compute, sample, vcf_file, "sum"])
    subprocess.call(["rm", vcf_file])

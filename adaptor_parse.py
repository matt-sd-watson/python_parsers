import os
from os import listdir
from os.path import isfile, join
from itertools import islice
import gzip
import pandas as pd
from pathlib import Path
import glob
import fnmatch


def adaptor_parse(filepath, output_file):
    filepaths = []
    fastqs = []
    for root, dirs, files in os.walk(os.path.abspath(filepath)):
        for file in files:
            filepaths.append(os.path.join(root, file))
            fastqs.append(file)

    filepaths = [name for name in filepaths if 'fastq.gz' in name]
    fastqs = [name for name in fastqs if 'fastq.gz' in name]

    adaptors = []

    for name in filepaths:
        with gzip.open(name, 'rt') as handle:
            for line in islice(handle, 0, 1):

                adaptors.append(line[58:66])

    cut_adapt = []
    for i, h, z in zip(adaptors, fastqs, filepaths):
            cut_adapt.append("cutadapt -a {} -o /Users/mattsdwatson/VIB_proj_3425/fastq_pp/{} {}".format(i, h, z))

    with open(output_file, "w") as handle:
        for element in cut_adapt:
            handle.write(element)
            handle.write("\n")


# example of using this parsing function
adaptor_parse("/Users/mattsdwatson/VIB_proj_3425/fastq/", "cutadapt.sh")

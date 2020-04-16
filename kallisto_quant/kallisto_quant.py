# This script generate the shell command script for bulk quantification of rna seq libraries
# using kallisto quant with the single read option

from os import listdir
from os.path import isfile, join


# full file path names or fastqs were generated using
# find $(pwd) -maxdepth 1 -type f -not -path '*/\.*' > files_fullpath.txt

# partial file names (without the complete directory path) for the fastqs were generated using
# ls $search_path > filename.txt

def kallisto_script(fullnames, partial_names, output_file):
    fullnames_read = [f for f in open(fullnames, "r")]

    partial_name_read = [f for f in open(partial_names, "r")]

    kallisto_commands = []
    for i, k in zip(fullnames_read, partial_name_read):
        # shell command line for kallisto rna seq quantification of single read libraries
        kallisto_commands.append("kallisto quant -i index -o output/{} --single -l 200 -s 20 {}".format(k, i))

    with open(output_file, "w") as handle:
        # include the bash shebang line that is preferred for portability
        handle.write("#!/usr/bin/env bash")
        handle.write("\n")
        for element in kallisto_commands:
            handle.write(element)
            handle.write("\n")

# sample execution of the file function


kallisto_script("/Users/mattsdwatson/VIB_proj_3425/files_fullpath.txt",
                "/Users/mattsdwatson/VIB_proj_3425/filename.txt",
                "kallisto_quant.sh")

import re

# input file corresponds to all of the coding sequences for chromosome 2L
# of Drosophila melanogaster
file = "/Users/mattsdwatson/Documents/KU Leuven/Courses/" \
       "Practical Computing for Bioinformatics/example_exam/" \
         "drosophila_sequence.txt"

readfile = open(file, "r")
# enable the parsers to read the file
read = readfile.read()

# regex to find any multi-line piece of text that contains only ATCG characters, which
# stops at the headers that contain multi-character inputs other than nucleotides
sequences = re.findall(r"((?:[ATGC]+\n)+)", read)

# regex to find all gene names (pattern starts with gene and ends with ])
genes = re.findall("gene=.*?\]", read)

names = []
for gene in genes:
    # split the gene regex results to keep only the gene name
    # any time that the characters in brackets are found, perform a split
    split = re.split("([=\]]+)", gene)
    # keep only the third entry in the split corresponding to the complete gene name
    names += split[2::3]

# recover arginine frequencies in each coding sequence
frequencies = []
for seq in sequences:
    # read the nucleotides in groups of 3 for each line in the coding sequences
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    match = []
    for text in codons:
        # regex for argigine codons in each group of 3 from the codons
        # any character within [] can be used at the designated position
        match += re.findall("AG[AG]|CG[CGAT]", text)
    # save the number of times the match occurs in each section to frequencies
    frequencies += [len(match)]

# recover specific helix turn helix motifs in coding sequences
# if the motif is found at an index, append the gene name with the same index
genes_with_motifs = []
for i in range(0, len(sequences)):
    # check if the specific helix motif is found in the sequences list at index i
    # helix motif is of a specific length with designated codons found at predicted sites
    if re.search("[ATGC]{12}GC[ATCG]{7}CT[ATCG]GG[ATCG]{16}GT[ATCG]{7}", sequences[i]):
        # if the motif is found at i, append the gene name from name at i to genes with
        # motifs
        genes_with_motifs += [names[i]]

# read motif genes as single lines
genes_with_motifs = '\n'.join(genes_with_motifs)
print(genes_with_motifs)

# write all genes to file to use in a BioMart search
with open("genes.txt", "w") as write:
    for item in genes_with_motifs:
        write.write(item)

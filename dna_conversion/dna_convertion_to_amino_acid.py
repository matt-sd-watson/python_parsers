import re


def convert_coding_file(file, output_file=None):
    open_file = open(file, "r").read()
    sequences = re.findall(r"((?:[ATGC]+\n)+)", open_file)

# establish a dictionary of the amino acids and their corresponding abbreviations
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

# get the codons read in groups of three for translation
    codons = []
    for frame in sequences:
        for i in range(0, len(frame), 3):
            codon = frame[i:i + 3]
            codons.append(codon)
            if codon in ["TAG", "TAA", "TGA"]:
                break

    # establish a list of the proteins
    # if the protein in the table is present in the sequene, append to the protein list
    proteins = []
    for i in codons:
        if i in table:
            proteins.append(table.get(i))
        else:
            # append an X for any position where the amino acid cannot be read
            proteins.append("X")

    # concatenate all proteins into a continuous amino acid sequence
    final_sequence = ''.join(str(elem) for elem in proteins)
    individual_sequences = final_sequence.split("_")

    # if the output file is selected, write the amino acid sequences to file
    if output_file:
        with open(output_file, "w") as handle:
            for line in individual_sequences:
                handle.write(line)
                handle.write("\n")
    else:

        return individual_sequences





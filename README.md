# python_parsers
Python scripts for parsing various data files for extraction of relevant biological information

parse_PDB accepts a Protein Data Bank file format (such as 5kk.pdb), and parses for the distribution of atomic coordinates between hydrophilic and hydrophobic amino acids. 

parse_chromosome accepts a .txt coding sequence or FASTA file (such as drosophila_2L.txt), and retrieves all coding sequences with the corresponding genes names using matched indices. Parsing of a specific helix turn helix motif in the nucleotide sequences returns a list of candidate genes containing this motif, which can be used in a BioMart search for gene ontology. 

adaptor_parse accepts a file folder containing any fastq files held in within the sub-directories of the main input directory. It parses these files to retrieve their multiplexing adaptors for creating bulk scripts to execute cutadapt for adaptor trimming prior to alignment. 

kallisto_quant accepts two .txt files for input: one containing a list of partial filenames for rna-seq libraries, and another containing the list of absolute/full paths for the same libraries. it generates a shell script that will quantify all of the libraries using kallisto quant in the single read format. 

dna_conversion accepts a nucleotide coding file and returns all open reading frames
with the amino acid equivalent codons. 

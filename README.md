# python_parsers
Python scripts for parsing various data files for extraction of relevant biological information

parse_PDB accepts a Protein Data Bank file format (such as 5kk.pdb), and parses for the distribution of atomic coordinates between hydrophilic and hydrophobic amino acids. 

parse_chromosome accepts a .txt coding sequence or FASTA file (such as drosophila_2L.txt), and retrieves all coding sequences with the corresponding genes names using matched indices. Parsing of a specific helix turn helix motif in the nucleotide sequences returns a list of candidate genes containing this motif, which can be used in a BioMart search for gene ontology. 

adaptor_parse accepts a file folder containing any fastq files held in the main directory of any sub-directory. it parses these files to retrieve their multiplexing adaptors for creating bulk scripts to execute cutadapt for adaptor trimming prior to alignment. 

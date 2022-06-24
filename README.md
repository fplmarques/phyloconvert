This script parses NCBI nucleotide data in GenBank full format
and output fasta file with different headings format.


usage: phyloconvert.py [-h] -i INPUT_FILE [-o OUTPUT_FORMAT] [-m MAFFT_ALIGNMENT] [-s STATS]

Convert files (GenBank or FAST) for phylogenetic analyses and of sequence statistics.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Input file name *.gb for GenBank full format or *.fas for FASTA format (default: None)
  -o OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
                        'n' for accession numbers, 't' for taxon names, and 'b' for both (default: b)
  -m MAFFT_ALIGNMENT, --mafft_alignment MAFFT_ALIGNMENT
                        Make mafft alignment (default: False)
  -s STATS, --stats STATS
                        Provide sequence statistics (default: False)

This script parses NCBI nucleotide data in GenBank full format
and output fasta file with different headings format.

## Help instuctions

usage: phyloconvert.py [-h] -i INPUT_FILE [-o OUTPUT_FORMAT] [-m MAFFT_ALIGNMENT] [-s STATS]

Convert files (GenBank or FAST) for phylogenetic analyses and of sequence statistics.

#### Optional arguments:
-h, --help            show this help message and exit

-i INPUT_FILE, --input_file INPUT_FILE
Input file name *.gb for GenBank full format or *.fas for FASTA format (default: None)

-o OUTPUT_FORMAT, --output_format OUTPUT_FORMAT
'n' for accession numbers, 't' for taxon names, and 'b' for both (default: b)

-m MAFFT_ALIGNMENT, --mafft_alignment MAFFT_ALIGNMENT
Make mafft alignment (default: False)

-s STATS, --stats STATS
Provide sequence statistics (default: False)


## Examples

### Ouputs unaligned sequences using NCBI accession numbers:
> ./phyloconvert -i test.gb -o n

Other options for -o are to incluce only taxon names (t) and the default option is to inlude both (b).

### Outputs aligned sequences:
> ./phyloconvert -i test_taxa.fas -m=True

This command line alignes **test_taxa.fas** (requires FASTA format). Alignment calls require [MAFFT](https://mafft.cbrc.jp/alignment/software/) to be installed in your system. MAFFT will use `--maxiterate 1000 --globalpair` to align sequences.

### Adding sequence basic statistics:
> ./phyloconvert -i test_taxa.fas -m=True -s=True

Same as above, but including a standart output with the following summary statistics:

| Parameter      | Value |
| ----------- | -----------: |
| Number of sequences: | 32 |
| Length of aligned sequences: | 1689 |
| Max. unaligned sequence: | 1553 |
| Min. unaligned sequence: | 975 |

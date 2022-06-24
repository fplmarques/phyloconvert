#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##
# LIBRARIES
##

import os
import sys
from Bio import SeqIO
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

##
# FUNCTIONS
##

def check_input_type(in_file):
    """
      --> This function check extension of the input file to define de behavious of the program
    :param in_file: input file called after option -i
    :return: extension type
    """
    input_type = in_file.split('.')
    global prefix_file_name
    prefix_file_name = input_type[0]
    if input_type[1] == 'gb':
        return 'gb'
    if input_type[1] == 'fa' or input_type[1] == 'fas' or input_type[1] == 'fasta':
        return 'fa'


def parse_gb_file(gb_file):
    """
      --> Function parses gb files
    :param gb_file: NCBI records in genbank full format (*.gb)
    :return: dictionary with accession numbers, taxon names and sequence data
    """
    input_file = open(gb_file)
    accession_numbers = list()
    names = list()
    sequence_data = list()
    global n_records
    n_records = 0
    for seq_record in SeqIO.parse(input_file, "genbank"):
        accession_numbers.append(seq_record.name)
        names.append(seq_record.annotations["organism"].replace(' ', '_').replace('.', ''))
        sequence_data.append(str(seq_record.seq))
        n_records += 1
    input_file.close()
    return {'ncbi_number': accession_numbers, 'taxon': names, 'seqs': sequence_data}


def parse_fasta_file(fasta_file):
    """
      --> Function parses FASTA files
    :param fasta_file: Any sequence file in FASTA format
    :return: Dictionay with sequence IDs and sequence data
    """
    input_file = open(fasta_file)
    sequence_id = list()
    sequence_data = list()
    global n_records
    n_records = 0
    for seq_record in SeqIO.parse(input_file, "fasta"):
        sequence_id.append(seq_record.id)
        sequence_data.append(str(seq_record.seq))
        n_records += 1
    input_file.close()
    return {'organism': sequence_id, 'seqs': sequence_data}


def print_fasta(dict, n_records, in_format='gb', out_format='b', tmp=True):
    """
      --> Print sequence file according to selected options
    :param dict: Dictionary provided by either parse_gb_file(gb_file) or parse_fasta_file(fasta_file)
    :param n_records: Nunber of sequences collected by functions that generate dictionaries
    :param in_format: Input format, either Genbank full (*.gb) or FASTA (*.fa|*.fas|*.fasta)
    :param out_format: Define formats for sequence IDs: n=accession numbers, t=taxa, and b=both
    :param tmp: Boolean for writing tmp.tmp
    :return: FASTA file
    """
    if in_format == 'gb':
        if tmp is True:
            with open('tmp.tmp', 'w') as tmpf:
                for count in range(0, n_records):
                    if out_format == 'n':
                        tmpf.write(f'>{dict["ncbi_number"][count]}\n{dict["seqs"][count]}\n')
                    elif out_format == 't':
                        tmpf.write(f'>{dict["taxon"][count]}\n{dict["seqs"][count]}\n')
                    else:
                        tmpf.write(f'>{dict["ncbi_number"][count]}_{dict["taxon"][count]}\n{dict["seqs"][count]}\n')
        else:
            for count in range(0, n_records):
                if out_format == 'n':
                    print(f'>{dict["ncbi_number"][count]}\n{dict["seqs"][count]}')
                elif out_format == 't':
                    print(f'>{dict["taxon"][count]}\n{dict["seqs"][count]}')
                else:
                    print(f'>{dict["ncbi_number"][count]}_{dict["taxon"][count]}\n{dict["seqs"][count]}')
    else:
        if tmp is True:
            with open('tmp.tmp', 'w') as tmpf:
                for count in range(0, n_records):
                    tmpf.write(f'>{dict["organism"][count]}\n{dict["seqs"][count]}\n')
        else:
            for count in range(0, n_records):
                print(f'>{dict["organism"][count]}\n{dict["seqs"][count]}')


def mafft_align(input_file, output_file):
    """
      --> Align sequences (FASTA format) using MAFFT
    :param input_file: name of FASTA file of unaligned sequences
    :param output_file: name of output FASTA file of aligned sequences
    :return: FASTA file of aligned sequences
    """
    os.system(f'mafft --maxiterate 1000 --globalpair ./{input_file} > {output_file}')


def table_head(msg,col='white',bkg='dark_gray'):
    """
      --> Print the head of table
    :param msg: Centered titled
    :param col: font color
    :param bkg: backgroup color
    :return:
    """
    color = {'dark_gray':90, 'red':91, 'green':92, 'yellow':93, 'white':98}
    bkgd = {'gray':47, 'white':48, 'dark_gray':49, 'red':31, 'green':32, 'yellow':33}
    global size
    size = len(msg)+30
    print('-'*size)
    msg = msg.center(size).upper()
    print(f'\033[{color[col]}:{bkgd[bkg]}m{msg}\033[00m')
    print('-'*size)


def sequence_stats(dict):
    """
      --> Get stattistics from aligned FASTA file
    :param dict: returned dictionaries either parse_gb_file(gb_file) or parse_fasta_file(fasta_file)
    :return: Print stats in table for which the head is printed by table_head()
    """
    n_records = len(dict['organism'])
    aligned_seq_len = len(dict['seqs'][0])
    min_len = max_len = len(dict['seqs'][0].replace('-', ''))
    for count in range(0, n_records):
        get_ungapped_len = len(dict['seqs'][count].replace('-', ''))
        if get_ungapped_len < min_len and count != 0:
            min_len = get_ungapped_len
        elif get_ungapped_len > max_len and count != 0:
            max_len = get_ungapped_len
    table_head('SUMMARY STATISTICS', col='green', bkg='dark_gray')
    print(f' Number of sequences: \t\t\t {n_records}\n'
          f' Length of aligned sequences: \t\t {aligned_seq_len}\n'
          f' Max. unaligned sequence: \t\t {max_len}\n'
          f' Min. unaligned sequence: \t\t {min_len}')
    print('-'*size)


##
# MAIN FUNCTION
##


def main(args):
    # Setting parameters
    in_file = args['input_file']
    out_format = args['output_format']
    mafft = args['mafft_alignment']
    stats = args['stats']
    if stats:
        mafft = True

    # main program
    file_type = check_input_type(in_file)
    if file_type == 'fa':
        dict = parse_fasta_file(in_file)
        print_fasta(dict, n_records, in_format='fa', tmp=True)
        print(f'I am not sure what you want me to do here! Use option "-m True" for alignment '
              f'and "-s True" for sequence statistics.')
        if mafft and not stats:
            mafft_align('tmp.tmp', f'{prefix_file_name}_aln.fas')
            os.system('rm -f tmp.tmp')
        if stats and mafft:
            if os.path.exists(f'./{prefix_file_name}_aln.fas'):
                aln = parse_fasta_file(f'./{prefix_file_name}_aln.fas')
                sequence_stats(aln)
            else:
                mafft_align('tmp.tmp', f'{prefix_file_name}_aln.fas')
                aln = parse_fasta_file(f'./{prefix_file_name}_aln.fas')
                sequence_stats(aln)
    elif file_type == 'gb':
        dict = parse_gb_file(in_file)
        if out_format == 'b':
            print_fasta(dict, n_records, in_format='gb', out_format='b', tmp=True)
            os.system(f'mv tmp.tmp {prefix_file_name}_acc_and_name.fas')
            if mafft and not stats:
                mafft_align(f'{prefix_file_name}_acc_and_name.fas', f'{prefix_file_name}_acc_and_name_aln.fas')
            if stats and mafft:
                if os.path.exists(f'./{prefix_file_name}_acc_and_name_aln.fas'):
                    aln = parse_fasta_file(f'./{prefix_file_name}_acc_and_name_aln.fas')
                    sequence_stats(aln)
                else:
                    mafft_align(f'{prefix_file_name}_acc_and_name.fas', f'{prefix_file_name}_acc_and_name_aln.fas')
                    aln = parse_fasta_file(f'./{prefix_file_name}_acc_and_name_aln.fas')
                    sequence_stats(aln)
        elif out_format == 'n':
            print_fasta(dict, n_records, in_format='gb', out_format='n', tmp=True)
            os.system(f'mv tmp.tmp {prefix_file_name}_acc_num.fas')
            if mafft and not stats:
                mafft_align(f'{prefix_file_name}_acc_num.fas', '{prefix_file_name}_acc_num_aln.fas')
            if stats and mafft:
                if os.path.exists(f'./{prefix_file_name}_acc_num_aln.fas'):
                    aln = parse_fasta_file(f'./{prefix_file_name}_acc_num_aln.fas')
                    sequence_stats(aln)
                else:
                    mafft_align(f'{prefix_file_name}_acc_num.fas', f'{prefix_file_name}_acc_num_aln.fas')
                    aln = parse_fasta_file(f'./{prefix_file_name}_acc_num_aln.fas')
                    sequence_stats(aln)
        elif out_format == 't':
            print_fasta(dict, n_records, in_format='gb', out_format='t', tmp=True)
            os.system(f'mv tmp.tmp {prefix_file_name}_taxa.fas')
            if mafft and not stats:
                mafft_align(f'{prefix_file_name}_taxa.fas', f'{prefix_file_name}_taxa_aln.fas')
            if stats and mafft:
                if os.path.exists(f'./{prefix_file_name}_taxa_aln.fas'):
                    aln = parse_fasta_file(f'./{prefix_file_name}_taxa_aln.fas')
                    sequence_stats(aln)
                else:
                    mafft_align(f'{prefix_file_name}_taxa.fas', f'{prefix_file_name}_taxa_aln.fas')
                    aln = parse_fasta_file(f'./{prefix_file_name}_taxa_aln.fas')
                    sequence_stats(aln)
    else:
        print(f'\033[91m ERROR: Input file format not supported! \033[00m')

##
# INITIALIZATION
##

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
                    description='Convert files (GenBank or FAST) for phylogenetic analyses and of sequence statistics.',
                    epilog=f'\033[92m Fernando P.L. Marques June, 2022 \033[00m')
    parser.add_argument('-i', '--input_file', type=str, required=True,
                    help='Input file name *.gb for GenBank full format or *.fas for FASTA format')
    parser.add_argument('-o', '--output_format', default='b', type=str,
                    help="'n' for accession numbers, 't' for taxon names, and  'b' for both")
    parser.add_argument('-m', '--mafft_alignment', default=False, type=bool,
                    help="Make mafft alignment")
    parser.add_argument('-s', '--stats', default=False, type=bool,
                    help='Provide sequence statistics')
    args = vars(parser.parse_args())
    main(args)  # This will call the main function

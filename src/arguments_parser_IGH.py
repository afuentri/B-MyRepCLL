# ! /usr/env/python3.5                      #
# Azahara Maria Fuentes Trillo              #
# Unidad de Genomica y Diagnostico Genetico #
# INCLIVA                                   #
# 2018-04-11                                #
#############################################

# Argsparse in standard mode for general pipeline
# Import modules
import argparse
import sys

def parse_args():

    """Parse args with argparse"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--pipeline', '-m', dest='pipeline', required=False, help='Run pipeline in mode analysis', action='store_true')
    parser.add_argument('--test', '-t', dest='test', required=False, help='Run pipeline in mode test with a clonal control sample.', action='store_true')
    parser.add_argument('--input', '-f', dest='inputf', required='--pipeline' in sys.argv, help='path to input FASTQ files folder')
    parser.add_argument('--exclude', '-e', dest='exclude', required=False, help='Directory to exclude from fastq seeking')
    parser.add_argument('--sanger', '-s', dest='sanger', required=False, help='Select this option when using sanger sequences', action='store_true')
    parser.add_argument('--clonality', '-c', dest='clonality', required=False, help='Select this option when exhaustive clonality analysis needs to be performed (slower, reporting results for each read)', action='store_true')
    parser.add_argument('--primers', '-a', dest='primers', required=False, help='Select this option to trim primers giving a FASTA file. For both sides trimming include the name "left" and "right" in your primers files, so as the name of the right primers can be inferred')
    parser.add_argument('--bothsides', '-i', dest='sides', required=False, help='Select this option to trim primers to the left and the right of the sequences. FASTA files should be given with the option --primers', action='store_true')
    parser.add_argument('--basal', '-b', dest='basal', required=False, help='Select this option to use a default basal filter of 3.6 percent or use as basal filter IGHV clone with the highest representation in a policlonal sample if samples containing "policlonal" in sample name are given. If this argument is not given there will be no filters applied.', action='store_true')
    parser.add_argument('-d', '--cdr3s', help='annotate cdr3s over 100 reads', action='store_true')
    parser.add_argument('--out', '-o', dest='dest', required=True, help='Destiny folder to store the results')
    parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('--processes', '-p', dest='processes', required='--pipeline' in sys.argv, help='Number of processes')
    parser.add_argument('--overwrite', '-w', dest='owrt', required=False, help='If given, all files excepting merged FASTQ files (raw data must be removed by the user if wanted) will be overwritten.', action='store_true') 
    
    #parser.add_argument('-c', '--clust', help='clustering output repository (if performed)', action="store_true")
    args = parser.parse_args()

    return args

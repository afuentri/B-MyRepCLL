# ! /usr/env/python                         #
# Azahara Maria Fuentes Trillo              #
# Unidad de Genomica y Diagnostico Genetico #
# INCLIVA                                   #
# 2017-04                                   #
#############################################

# Comments here

"""This script performs the generation of two tables in .csv format. 
The first one is probable_(regionname)regions.csv, with the counts of references 
with more than one read mapped against it for all samples. 
The other table called "resume_(regionname)regions.csv", with the 10 references 
per sample which have the maximum number of reads mapped.
"""

# Import modules
import argparse
import string
import os
import logging as log
import heapq
import subprocess

def arguments_parser():

    """Parse args with argparse"""

    parser = argparse.ArgumentParser(description='IGH probable regions', argument_default=None)
    parser.add_argument('--fof_stats', '-f', dest='FofStats', required=True, help='Fof for samtools idx stats output')
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument('--out', '-o', dest='f_out', required=True, help='File for probable regions output')
    
    args = parser.parse_args()

    return args


def verbose(fichero_out):

    """Verbosity"""

    log.info('Generating file %s', fichero_out)


def read_fof(input):

    """Open file and store as text"""
    stats_files = []
    while True:
        try:
            f_hand = open(input, 'r')
        except:
            raise IOError(input + ' is not a file. Try again')
            continue
        else:
            log.info('Reading %s succesfully', input)
            # Go through fof
            for line in f_hand:
                line = line.strip()
                stats_files.append(line)

            f_hand.close()
            break

    return stats_files


def extract_fragmentreads(ref, bam):

    """"""
    dfrag = {'leader': [], 'FR1': [], 'FR2': [], 'FR3': []}
    CMD = 'samtools view {} -F4 {}'.format(bam, ref)
    proce = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE)
    out = proce.communicate()[0].decode("utf-8").split('\n')
    
    for line in out:
        
        if line != '':
            line = line.strip()
            
            seqid = line.split()[0]
            pos = int(line.split()[3])
            
            if pos < 120 and pos >= 40:
                dfrag['FR1'].append(seqid)
            elif pos < 200 and pos >= 120:
                dfrag['FR2'].append(seqid)
            elif pos < 280 and pos >= 200:
                dfrag['FR3'].append(seqid)
            elif pos < 40 and pos >= 0:
                dfrag['leader'].append(seqid)
        
    return dfrag


def get_regions(ini, name, f_reg, f_res, bam):

    """If a region has more or equal to one read mapped to it, it is reported"""

    reg = {}
    for region in ini:
        
        region = region.strip()
        if region.split('\t')[0] != '*':
            ref = region.split("\t")[0].replace('/','-')
            ref_len = region.split("\t")[1]
            reads_mapped = int(region.split("\t")[2])
            reg[ref] = reads_mapped
            dfrag = extract_fragmentreads(ref, bam)
        
            if reads_mapped > 0:

                f_reg.write(str(name) + "," +  str(ref) + "," + str(reads_mapped) + "," + str(ref_len)
                            + "," + str(len(dfrag['leader'])) + "," + str(len(dfrag['FR1'])) + ","
                            + str(len(dfrag['FR2'])) + "," + str(len(dfrag['FR3'])) +'\n')

    ## Extract the first 10 maximums
    A = heapq.nlargest(10, reg, key=reg.get)

    for region in A:
        if reg[region] > 0:
            dfrag2 = extract_fragmentreads(region, bam)
            f_res.write(str(name) + "," + str(region) + "," + str(reg[region]) +
                        "," + str(len(dfrag2['leader'])) + "," + str(len(dfrag2['FR1'])) + "," +
                        str(len(dfrag2['FR2'])) + "," + str(len(dfrag2['FR3'])) + "\n")
        
        
    
def main():

    """The main program"""

    # Get arguments
    args = arguments_parser()

    #   Logging and verbose
    if args.verbose:
        log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)
    else:
        log.basicConfig(format='%(levelname)s:%(message)s', level=log.WARNING)

    # Print header
    f_regions = open(args.f_out, 'w')
    f_regions.write("sample_name,region,n_reads_mapped,region_length,reads_leader,reads_FR1,reads_FR2,reads_FR3\n")

    ## On resume file as well
    res = args.f_out.replace('probable', 'resume')
    f_resume = open(res, 'w')

    f_resume.write("sample_name,region,n_reads_mapped,reads_leader,reads_FR1,reads_FR2,reads_FR3\n")

    # Open fof
    stats_list = read_fof(args.FofStats)

    # Open idxstat files
    for sample in stats_list:

        bam = sample.replace('statsV', 'bamsV').replace('statsJ', 'bamsJ').replace('-stats.txt', '-sorted.bam')
        stats = open(sample, 'r')
        # Sample name
        sample_name = os.path.basename(sample).replace("-stats.txt", "")
        log.info('Parsing %s sample stat file', sample_name)

        # Select regions with more or equal than one read mapped to an allele and the resume with the maximums
        get_regions(stats, sample_name, f_regions, f_resume, bam)

        stats.close()

    # Close files

    f_regions.close()
    f_resume.close()

if __name__ == "__main__":
    main()

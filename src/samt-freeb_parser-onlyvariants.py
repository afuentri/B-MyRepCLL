# ! /usr/env/python                         #
# Azahara Maria Fuentes Trillo              # 
# Veronica Lendinez Tortajada               #
# Unidad de Genomica y Diagnostico Genetico #
# INCLIVA                                   #
# 2018-06-01                                #
#############################################

# Comments here
#
# Import modules
import argparse
import os
import logging as log
import string
import gzip
import ast
import sys

def arguments_parser():

    """Parse args with argparse"""

    parser = argparse.ArgumentParser(description='VCF parser', argument_default=None)
    parser.add_argument('--fof', '-V', dest='FofVcf', required=True, help='Fof for whole genome VCFs')
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

    args = parser.parse_args()

    return args


def extract_vcf_list(input):

    """"""
    # Open FOF
    while True:
        try:
            fof = open(input, 'r')
        except:
            raise IOError(input + ' is not a file. Try again')
            continue 
        else:
            log.info('Reading %s succesfully', input)
            
            list = fof.readlines()
            list = map(str.strip, list)
            
            return list
            
            break


def sample_info(list):

    """"""
    log.info('Samples specified: ')
    
    for sample in list:
        log.info('%s', os.path.basename(sample))
    

def parse_vcf(f_hand, sample):

    """"""
    vcf_file = gzip.open(f_hand, 'r')

    for line in vcf_file:
        if not line.startswith('#'):
    
            line = line.strip().split('\t')
            chr, position, ref, alt = (line[0], line[1], line[3],
                                       line[4].replace(',',';'))
                    
            if alt != '.':
                #rs = line.split('\t')[2]
                form = dict(zip(line[8].split(':'), line[9].split(':')))
                GT = form['GT'].replace('\n', '')
                DP = form['DP'].replace('\n', '')
                AD = form['AD'].replace(',',';').replace('\n', '')
                ## VAF
                VAF = []
                for k in range(0, len(alt.split(','))):
                               
                    VAF.append(str(round(float(AD.split(';')[k + 1])/float(DP), 2)))
                               
                if 'ADF' in form:
                    ADF = form['ADF'].replace(',',';').replace('\n', '')
                else:
                    ADF = ''
                if 'ADR' in form:
                    ADR = form['ADR'].replace(',',';').replace('\n', '')
                else:
                    ADR = ''
                if 'GQ' in form: 
                    GQ = form['GQ'].replace('\n', '')
                else:
                    GQ = ''
                if 'RO' in form:
                    RO = form['RO'].replace(',',';').replace('\n', '')
                else:
                    RO = ''
                if 'AO' in form:
                    AO = form['AO'].replace(',',';').replace('\n', '')
                else:
                    AO = ''

                print(str(sample.replace('.vcf.gz', '')) + ',' + str(chr) + ',' + str(position) +\
                      ',' + str(ref) + ',' + str(alt) + ',' + str(GT) + ',' + str(DP) + \
                      ',' + str(AD) + ',' + ';'.join(VAF) + ',' + str(ADF) +\
                      ',' + str(ADR) + ','+ str(GQ) + ',' + str(RO) + ',' + str(AO))

    vcf_file.close()
                      

def main():

    
    """The main program"""
    
    # Get arguments
    args = arguments_parser()
    # Logging and verbose
    if args.verbose:
        log.basicConfig(format='%(levelname)s: %(message)s', level=log.DEBUG)
    else:
        log.basicConfig(format='%(levelname)s: %(message)s', level=log.WARNING)
        
    log.info('Starting VCF parsing script ')
    
    # Read fof
    VcfList = extract_vcf_list(args.FofVcf)
    
    # Print sample info
    sample_info(VcfList)

    # Parse VCF
    print("sample,Vallele,position,ref,alt,GT,DP,AD,VAF,ADF,ADR,GQ,RO,AO")
    
    
    for vcf in VcfList:
        parse_vcf(vcf, os.path.basename(vcf))

        
if __name__ == "__main__":
    main()

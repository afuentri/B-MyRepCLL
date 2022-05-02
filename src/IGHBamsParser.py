# ! /usr/env/python                         #
# Azahara Maria Fuentes Trillo              #
# Unidad de Genomica y Diagnostico Genetico #
# INCLIVA                                   #
# init:2018-04                              #
#############################################

# Comments here
"""This script parses V, D, and  sorted bams for each sample given in a FOF and 
annotates each read of V BAM with the mapping information from the other two bams
to find a complete subtype for all regions
"""

# Import modules
import argparse
import os
import logging as log
import pysam
import string

def arguments_parser():

    """Parse args with argparse"""

    parser = argparse.ArgumentParser(description='IGH mapping info', argument_default=None)
    parser.add_argument('--Vfof', '-V', dest='FofBamsV', required=True, help='Fof for reads mapped against V region bams')
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument('--out', '-o', dest='out', required=True, help='Output directory')
    parser.add_argument("-c", "--clonal", help="run clonality mode", action="store_true")
    args = parser.parse_args()

    return args


def verbose(fichero_out):

    """Verbosity"""

    print('Generating file %s') %fichero_out


def read_bam(bam):

    """Open bam file"""

    while True:
        try:
            bamfile = pysam.AlignmentFile(bam, "rb")
        except:
            raise IOError(bam + ' is not a file. Try again')
            continue
        else:
            log.info('Reading %s succesfully', bam)

        return bamfile


    return bamfile


def read_fof(input):

    """Open file and store as text"""
    bam_files = []
    while True:
        try:
            f_hand = open(input, 'r')
        except:
            raise IOError(input + ' is not a file. Try again')
            continue #probar
        else:
            log.info('Reading %s succesfully', input)
            # Go through fof
            for line in f_hand:
                line = line.strip()
                bam_files.append(line)

            f_hand.close()
            break

    return bam_files


def parse_bam_file(bam):

    """Extract concrete information from a BAM file"""
    d = {}
    for read in bam.fetch(until_eof=True):
        # fetch() will only iterate over alignments in the SAM/BAM file so we use until_eof
        # when check_sq is true it defines chromosome names
        read_id = read.qname
        if read.tid != -1:
            region = bam.get_reference_name(read.tid)
            
        else:
            region = 'None'
    
        cigar = read.cigarstring
        sequence = read.seq
        
        if read_id not in d:
            d[read_id] = [sequence, cigar, region]
        else:
            d[read_id][0] = str(d[read_id][0]) + ';' + str(sequence)
            d[read_id][1] = str(d[read_id][1]) + ';' + str(cigar)
            d[read_id][2] = str(d[read_id][2]) + ';' + str(region)
            
    return d


def write_to_file(sample, f, info_V, info_J, *args, **kwargs):

    """Write info to file on TXT format"""
    # clonality files
    fv = kwargs.get('Vfile', None)
    fj = kwargs.get('Jfile', None)
    consensus_out = kwargs.get('cout', None)
    for sequence in info_V:
        v1 = None
        j1 = None
        if fv != None and fj != None:
            if 'None' not in info_V[sequence][2]:
                if len(info_V[sequence][2].split(';')) > 1:
                    
                    a1, a2 = info_V[sequence][2].split(';')[0:2]
                    if a1 == a2:
                        fv.write('%s_%s,%s,1,\n' %(sample, sequence, a1))
                        v1 = a1
                        v2 = None
                    else:
                        fv.write('%s_%s,%s,1,\n' %(sample, sequence, a1))
                        fv.write('%s_%s,%s,1,\n' %(sample, sequence, a2))
                        v1 = a1
                        v2 = a2
                        
                else:
                    fv.write('%s_%s,%s,1,\n' %(sample, sequence, info_V[sequence][2]))
                    v1 = info_V[sequence][2]
                    v2 = None
                    
            if 'None' not in info_J[sequence][2]:
                if len(info_J[sequence][2].split(';')) > 1:
                    b1, b2 = info_J[sequence][2].split(';')
                    if b1 == b2:
                        fj.write('%s_%s,%s,1\n' %(sample, sequence, b1))
                        j1 = b1
                        j2 = None
                    else:
                        fj.write('%s_%s,%s,1\n' %(sample, sequence, b1))
                        fj.write('%s_%s,%s,1\n' %(sample, sequence, b2))
                        j1 = b1
                        j2 = b2
                        
                else:
                    fj.write('%s_%s,%s,1\n' %(sample, sequence, info_J[sequence][2]))
                    j1 = info_J[sequence][2]
                    j2 = None
                    
            ## annotate        
            if v1 and j1:
                elementsv = [v1, v2]
                elementsj = [j1, j2]
                    
                elev = [x for x in elementsv if x]
                elej = [x for x in elementsj if x]
                    
                for e in elev:
                    for w in elej:
                        cname = '{}_{}_{}_{}-fb.fa'.format(sample, sequence.replace('-', '123456789').replace(':', '-'), e, w)
                        consensus_path = os.path.join(consensus_out, cname)
                        fconsensus = open(consensus_path, 'w')
                        fconsensus.write('>{}_{}\n{}\n'.format(e, w, info_V[sequence][0]))
                        fconsensus.close()
                    
        # info_files
        f.write('[ID] %s \n' %sequence)
        f.write('[sequence]\nV\t%s\nD\t%s\nJ\t%s\n' %(info_V[sequence][0], \
        None, info_J[sequence][0]))
        f.write('[feature]\tV\tD\tJ\n')
        f.write('[cigar]\t%s\t%s\t%s\n' %(info_V[sequence][1], \
        None, info_J[sequence][1]))
        f.write('[region]\t%s\t%s\t%s\n' %(info_V[sequence][2], \
        None, info_J[sequence][2]))
        f.write("-" * 100 + '\n')
        

def main():

    """The main program"""

    # Get arguments
    args = arguments_parser()

    # Logging and verbose
    if args.verbose:
        log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)
    else:
        log.basicConfig(format='%(levelname)s:%(message)s', level=log.WARNING)

    # Open files
    bam_list = read_fof(args.FofBamsV)
    samples = [os.path.basename(bam_list[i].replace('-sorted.bam','')) for i in range(len(bam_list))]
    log.info('Samples specified: %s', samples)
    
    # Set out dir
    out_dir = args.out

    # argument clonal
    if args.clonal:
        fresults = os.path.dirname(os.path.dirname(out_dir))
        consensus_out = os.path.join(fresults, 'consensus_complete')
        
        pv = os.path.join(out_dir, 'regions_clonalityV.csv')
        pj = os.path.join(out_dir, 'regions_clonalityJ.csv')
        pv_out = open(pv, 'w')
        pj_out = open(pj, 'w')
        pv_out.write('sample_name,region_max_reads_mapped,max_reads_mapped,length\n')
        pj_out.write('sample_name,region_max_reads_mapped,max_reads_mapped\n')        
        
    # Open bams
    for bam in bam_list:
        
        # Open output files
        sample_name = os.path.basename(bam.replace('-sorted.bam',''))
        path_bam = '/'.join(bam.split('/')[:-2])
        f = os.path.join(out_dir, 'info_bams_' + sample_name + '.txt')
        f_out = open(f, 'w')
        log.info('Parsing %s sample BAM files', sample_name)
        bamV = read_bam(bam)
        bamJ = read_bam(os.path.join(path_bam, 'bamsJ', os.path.basename(bam)))
        
        # Parse bam files
        info_V = parse_bam_file(bamV)
        info_J = parse_bam_file(bamJ)
       
        if args.clonal:
            write_to_file(sample_name, f_out, info_V, info_J,
                          cout=consensus_out, Vfile=pv_out,
                          Jfile=pj_out)
            
        else:
            write_to_file(sample_name, f_out, info_V, info_J)
                                  
        bamV.close()
        bamJ.close()

        f_out.close()
        
    if args.clonal:
        pv_out.close()
        pj_out.close()
        
if __name__ == "__main__":
    main()

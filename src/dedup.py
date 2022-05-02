#! /usr/bin/python
import gzip
import sys
fastq_file = sys.argv[1]

def read_file_simply_gzip(input):

    """Read file and return a list whith each line as element"""
    # Should remove empty lines here or in the main programme?
    while True:
        try:
            f_hand = gzip.open(input, 'r')
        except:
           raise IOError(input + ' is not a file. Try again')
           continue
        else:
            with gzip.open(input, 'rt', encoding='utf-8') as f_hand:
               list = f_hand.read().splitlines()

               f_hand.close()
               break
                  
    return list
                   


def uniq_fastq(fastq):

    """Remove duplicated reads from FASTQ file"""
    reads = []
    d1 = {}
    d2 = {}
    d_single = {}
    m = read_file_simply_gzip(fastq)
    
    # we use '|' because it is not a symbol encoding q-score
    for read in range(0, len(m), 4):
        id_read = m[read].split('/')[0]
        if len(m[read].split('/')) > 1:
            if int(m[read].split('/')[1]) == '1':
                if not id_read in d1:
                    d1[id_read] = '|'.join(m[read : read + 4])
            else:
                if not id_read in d2:
                    d2[id_read] = '|'.join(m[read : read + 4])

        else:
            if not id_read in d_single:
                d_single[id_read] = '|'.join(m[read : read + 4])
            
                    
    new = gzip.open(fastq, 'wt', encoding='utf-8')
    if d1 != {}:
        for read in d1:
            if read in d2:
                new.write(d1[read].replace('|','\n') + '\n')
                new.write(d2[read].replace('|','\n') + '\n')
        for read in d1:
            if not read in d2:
                new.write(d1[read].replace('|','\n') + '\n')
    if d2 != {}:
        for read in d2:
            if not read in d1:
                new.write(d2[read].replace('|','\n') + '\n')
    if d_single != {}:
        for read in d_single:
            new.write(d_single[read].replace('|','\n') + '\n')
            
    new.write('\n')
    new.close()

uniq_fastq(fastq_file)

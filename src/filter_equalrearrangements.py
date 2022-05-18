#! /usr/bin/python3

import sys
import os
import glob
import subprocess

## homology_table.csv
hom = sys.argv[1]
path_cons = sys.argv[2]

def execute(CMD):
    
    """"""
    subprocess.call([CMD], shell=True, executable='/bin/bash')


def read_file_simply(input):
     
    """Read file and return a list whith each line as element"""
    # Should remove empty lines here or in the main programme?
    while True:
        try:
            f_hand = open(input, 'r')
        except:
            raise IOError(input + ' is not a file. Try again')
            continue
        else:
            
            with open(input, encoding='utf-8') as f_hand:
                list = f_hand.read().splitlines()
       
            f_hand.close()
            break
       
    return list
         
    
def compare_consensus_equal(f1, f2):

    """"""
    n1 = [i for i in glob.glob(f1) if not 'FR3' in i]
    n2 = [i for i in glob.glob(f2) if not 'FR3' in i]
    with open(n1[0], 'r') as s1:
        seq1 = ''.join(s1.read().splitlines()[1:])
    with open(n2[0], 'r') as s2:
        seq2 = ''.join(s2.read().splitlines()[1:])    
    if seq1 == seq2:

        return True


def parse_emboss(fhand):

    """"""
    m = read_file_simply(fhand)
    if len(m) > 1:
        for e in range(len(m)):
            if m[e].startswith('# Length'):
                h = m[e].replace("#","") + ';' + m[e + 1].replace("#","")
                hom = h.split('(')[1].replace(")","").replace("%","")
            elif m[e].startswith('# Score'):
                score = m[e].replace("# Score: ","")
    else:
        h = ''
        hom = ''
        score = ''

    return hom, h, score


def compare_consensus(f1, f2):

    """"""
    n1 = [i for i in glob.glob(f1) if not 'FR3' in i]
    n2 = [i for i in glob.glob(f2) if not 'FR3' in i]
    print(f1, f2)
    if (n1 != []) and (n2 != []):
        
        out_path = 'temp.fa'
        CMD = ('water -asequence {} -bsequence {} -gapopen 10.0 '
               '-gapextend 0.5 -outfile {} '
               '-aformat pair').format(n1[0], n2[0], out_path)
        execute(CMD)

        identity, whole, score = parse_emboss(out_path)
        print(identity)
        if float(identity) >= 95:
            print(identity)
            return True    
    
    
do = hom.replace('.csv', 'simpalleles.csv')
dt = hom.replace('.csv', 'prepalleles.csv')
dtemp = open(dt, 'w')

d = {}
dall = {}
with open(hom, 'r') as fhand:

    for line in fhand:
        if not line.startswith('sample'):
            line = line.strip()
            sample = line.split(',')[0]
            IGHV = line.split(',')[1]
            nreads = line.split(',')[2]
        
            if not sample in d:

                d[sample] = {}
                d[sample][IGHV] = int(nreads)
                dall[sample] = {}
                dall[sample][IGHV] = line
            else:
                d[sample][IGHV] = int(nreads)
                dall[sample][IGHV] = line

            dtemp.write(','.join(line.split(',')[:2]) +
                        ',' + IGHV + ','  +
                        ','.join(line.split(',')[2:]) + '\n')
                
        else:
            dtemp.write('sample_name,Vregion,combined_alleles,reads_mapped,'
                        'reads_mapped_leader,reads_mapped_fr1'
                        ',reads_mapped_fr2,reads_mapped_fr3,'
                        'region_length,homology-IGHV_noFR3,'
                        'mutational_status_noFR3,joined,Jregion,'
                        'J_assigned,J_coincidence,IGHV-J,'
                        'consensus_length,IGHD_emboss,'
                        'insertions,deletions,ORF disruption'
                        ',majorproductive_seq'
                        ',nreads_majorproductive_seq,major_CDR3,'
                        'major_IGHD,major_prod\n')

print(d)
drem = {}
drep = {}
for e in d:

    print(e)

    ## sum total reads sample
    total_reads = sum(d[e].values())
    lalleles = list(d[e].keys())

    while len(lalleles) > 0:
        j = lalleles[0]
        #jgroup = j.split('-')[0]
        jgroup = '-'.join(j.split('-')[:2]).replace('D', '')
        
        lalleles.remove(j)
        
        for g in lalleles:

            if j != g:
                #ggroup = g.split('-')[0]
                ggroup = '-'.join(g.split('-')[:2]).replace('D', '')
                if ggroup == jgroup:
                    #diff = abs(d[e][j] - d[e][g])/float(total_reads)*100
                    #if diff < 0.7:
                    
                    name1 = '{}_{}*-fb.fa'.format(e, j)
                    name2 = '{}_{}*-fb.fa'.format(e, g)
                    cons1 = os.path.join(path_cons, name1)
                    cons2 = os.path.join(path_cons, name2)
                        
                    if compare_consensus(cons1, cons2):
                        if d[e][j] > d[e][g]:

                            maxallele = j
                            minallele = g
                        else:
                            maxallele = g
                            minallele = j
                            
                        print(j, g)#, diff)
                        if not e in drem:
                            drem[e] = {}
                        if not e in drep:
                            drep[e] = {}
                        if not minallele in drem:
                            drem[e][minallele] = maxallele
                        else:
                            if drem[e][minallele] in drem:
                                maxxallele = drem[e][minallele]
                                ## keep maxallele with bigger n reads
                                if int(dall[e][maxallele].split(',')[2]) > int(dall[e][maxxallele].split(',')[2]):
                                    drem[e][minallele] = maxallele
                                        
        
                                
                        if not maxallele in drep[e]:
                            drep[e][maxallele] = {}
                            drep[e][maxallele]['nreads'] = int(dall[e][j].split(',')[2]) + int(dall[e][g].split(',')[2])
                            drep[e][maxallele]['nlead'] = int(dall[e][j].split(',')[3]) + int(dall[e][g].split(',')[3])
                            drep[e][maxallele]['nfr1'] = int(dall[e][j].split(',')[4]) + int(dall[e][g].split(',')[4])
                            drep[e][maxallele]['nfr2'] = int(dall[e][j].split(',')[5]) + int(dall[e][g].split(',')[5])
                            drep[e][maxallele]['nfr3'] = int(dall[e][j].split(',')[6]) + int(dall[e][g].split(',')[6])
                            drep[e][maxallele]['comb'] = [maxallele, minallele]
                        else:
                            drep[e][maxallele]['nreads'] = int(drep[e][maxallele]['nreads']) + int(dall[e][minallele].split(',')[2])
                            drep[e][maxallele]['nlead'] = int(drep[e][maxallele]['nlead']) + int(dall[e][minallele].split(',')[3])
                            drep[e][maxallele]['nfr1'] = int(drep[e][maxallele]['nfr1']) + int(dall[e][minallele].split(',')[4])
                            drep[e][maxallele]['nfr2'] = int(drep[e][maxallele]['nfr2']) + int(dall[e][minallele].split(',')[5])
                            drep[e][maxallele]['nfr3'] = int(drep[e][maxallele]['nfr3']) + int(dall[e][minallele].split(',')[6])
                            drep[e][maxallele]['comb'].append(minallele)               
                            
dtemp.close()

dout = open(do, 'w')
with open(dt, 'r') as fhanr:
    
    for line in fhanr:
        if not line.startswith('sample'):
            line = line.strip()

            s = line.split(',')[0]
            v = line.split(',')[1]
            
            if s in drem:
                
                if v in drem[s]:
                    pass
                else:
                    if v in drep[s]:
                        l = '{},{},{},{},{},{},{},{}\n'.format(','.join(line.split(',')[:2]),
                                                               ';'.join(drep[s][v]['comb']), drep[s][v]['nreads'],
                                                               drep[s][v]['nlead'], drep[s][v]['nfr1'],
                                                               drep[s][v]['nfr2'], drep[s][v]['nfr3'],
                                                               ','.join(line.split(',')[8:]))
                        dout.write(l)
                    else:
                        dout.write(line + '\n')
            else:
                dout.write(line + '\n')
        else:
            dout.write(line)

dout.close()


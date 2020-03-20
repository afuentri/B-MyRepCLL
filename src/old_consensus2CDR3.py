#! /usr/bin/python3.5

"""This set of functions work to extract CDR3 in an intelligent way from the
consensus sequence"""

# import modules
import sys
import string
import re
import os
import pipeline
repo = os.environ['scripts_repo']
path_fastq_merge = os.path.join(repo, 'NGS-tools')
sys.path.append(path_fastq_merge)
import fastq_merge
from dir_management import *


def aminoacid_dictionary():

    """Codon dictionary to convert DNA sequence to protein"""
    d = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
         'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
         'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
         'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
         'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
         'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
         'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
         'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
         'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
         'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
         'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
         'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
         'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
         'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
         'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
         'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    return d


def read_FASTA(f):

    """"""
    with open(f) as fhand:
        return fhand.read().split('>')[1:]


def read_consensus(f):

    """Read the FASTA format consensus file"""
    return(seq.partition('\n') for seq in read_FASTA(f))


def read_sequence(f):

    """"""
    return [(info[1:], seq.replace('\n','').upper())
            for info, ignore, seq in read_consensus(f)]


def translate_codon(codon):

    """"""
    d = aminoacid_dictionary()
    return d.get(codon)


def drawback3(s, pos):

    """"""
    # generator
    aminoacid_p = ((translate_codon(s[p-3:p]), p) for p in range(pos, 3, -3))
    return aminoacid_p


def seek_sequence(pat, seq):

    return [m.start(0) for m in re.finditer(pat, seq)]


def sequence_remove(variants, seq):

    """"""
    new_seq = seq
    for i in variants: 
        pos, ref, alt = i
        seq.replace(alt, ref)
        post_pos = int(pos) + len(alt)
        pre_string = seq[:int(pos)-1]
        post_string = seq[post_pos-1:]
        new_seq = pre_string + ref + post_string
        
    return new_seq


def gather_variants(list_variants, sequence):

    """"""
    list_things = []
    for i in list_variants:
        alleles, pos = i.split(':')
        ref, alt = alleles.split('/')
        list_things.append((pos, ref, alt))

        sequence = sequence_remove(list_things, sequence)
    
    return sequence


def CDR3_minning(motif, min_cys, seq):

    """"""
    CDR3, CDR3_end, CDR3_start = ('', '', '')
    for m in motif:
        
        seq_CDR3 = []
        cys_codon = 1000
        aa = drawback3(seq, m + 3)
    
        while cys_codon > min_cys:
        
            cys_codon_pos = next(aa, None)
            if cys_codon_pos:
                ## 
                codon, cys_codon = cys_codon_pos
                seq_CDR3.append(codon)
                if codon == 'C':
                    break

        if CDR3 == '':
            CDR3, CDR3_end, CDR3_start = (None, None, None)
        else:
            CDR3 = ''.join(seq_CDR3[::-1])
            CDR3_end = m + 3
            CDR3_start = cys_codon -3
                    
                
    return CDR3, CDR3_start, CDR3_end
    #print seq_CDR3
    #CDR3 = ''.join(seq_CDR3[::-1])
    #CDR3_end = m
    #CDR3_start = cys_codon
    #return CDR3, (CDR3_start - 3), (CDR3_end + 3)


def get_CDR3(seq, mincys=80):

    """"""
    productive = False
    motif = seek_sequence('TGGGG[ACGT][ACGT]*?GG[ACGT]', seq)
    cys = seek_sequence('TG[TC]', seq)
    ## we have to check that there are sequences with cys
    ## greater than mincys
    ## we are removing all morif starting positions below 10
    list_cys = [ e for e in cys if e > mincys ]
    min_cys_big = [ e for e in list_cys if e > 10]
    min_motif = [ e for e in motif if e > 10]
    if min_motif == []:
    
        CDR3, start, end = (None, None, None)
    
    else:
        if min_cys_big != []:
            ## 80 is not valid for sequences in the middle
            min_cys = min(min_cys_big)
            CDR3, start, end = CDR3_minning(min_motif, min_cys, seq)
            productive = True
        else:
            CDR3, start, end = (None, None, None)
        
    
    return CDR3, start, end, productive


def decide_CDR3(A, B, C, P):
    
    """"""
    D_set = set(A)
    if not A[0] or not A[1]:

        return None, None, None, None
    
    elif len(D_set) > 1:

        if ('_' in A[0]) or (not A[0].endswith('W')):
            D = 1
               
        elif '_' in A[1] or not A[1].endswith('W'):
            D = 0
        else:
            D = A.index(max(A, key = len))
    else:
        D = 0
    
    return A[D], B[D], C[D], P[D]


def prepare_IGHD_blast(fasta, out):

    """"""
    out_blast = fasta.replace('.fasta', '.blast')
    pipeline.execute('blastn -db IMGT_D.fa -task blastn-short -evalue 1000 -dust no -strand plus'
                     ' -word_size 5 -query {} -outfmt 6 |  sort -k12 -nr > {}').format(fasta, out_blast)

    return out_blast


def prepare_IGHD_emboss(fasta, out, refD):

    """"""
    EW_list = []
    allelesD = pipeline.Vref_dictionary(refD)
    EW = os.path.join(out, os.path.basename(fasta).replace('.fasta', '-EW'))
    
    for k in allelesD:
        nameD = k.split('|')[1]
        out_file = '{}_{}.water'.format(EW, nameD)
        pipeline.execute('emboss_water -asequence {} -bsequence <(echo -e ">{}\\n{}"'
                         '-gapopen 10.0 -gapextend 0.5 -outfile {} -aformat pair').format(fasta, nameD,
                                                                                          allelesD[k], out_file)
        EW_list.append(out_file)
    ## parse alignment
    pipeline.parse_alignment_CDR3_emboss([fasta], out, EW_list)

    return EW_list, EW
    
    

def emboss_annotation(EW):

    """"""
    fblast = read_file_simply(EW)

    if fblast != []:
        if len(fblast) >= 3:
            CDR3 = ' '.join(fblast[0:3]).replace(',',' ')
        elif len(fblast) >=2:
            CDR3 = ' '.join(fblast[0:2]).replace(',',' ')
        elif len(fblast) >= 1:
            CDR3 = fblast[0].replace(',',' ')
    else:
        CDR3 = ''

    return CDR3


def blast_annotation(blast):
    
    """"""
    fblast = read_file_simply(blast)
    if fblast != []:
        if len(fblast) >= 3:
            CDR3 = (fblast[0].split('\t')[1].split('|')[1] + ' ' + ' '.join(fblast[0].split('\t')[2:4]) +
                    ';' + fblast[1].split('\t')[1].split('|')[1] + ' ' + ' '.join(fblast[1].split('\t')[2:4]) +
                    ';' + fblast[2].split('\t')[1].split('|')[1] + ' ' + ' '.join(fblast[2].split('\t')[2:4]))
        elif len(fblast) >=2:
            CDR3 = (fblast[0].split('\t')[1].split('|')[1] + ' ' + ' '.join(fblast[0].split('\t')[2:4]) +
                        ';' + fblast[1].split('\t')[1].split('|')[1] + ' ' + ' '.join(fblast[1].split('\t')[2:4]))
        elif len(fblast) >= 1:
            CDR3 = (fblast[0].split('\t')[1].split('|')[1] + ' ' + ' '.join(fblast[0].split('\t')[2:4]))
    else:
        CDR3 = ''

    return CDR3


def cdr3_extraction(fasta, ref_seq, VCF, disruption,
                    list_insertions, list_deletions, mincys):
    
    """
    SCRIPT FOR DETERMINING CDR3 IN A GIVEN SEQUENCE
    """
    
    ## main
    productivity = ''

    # revisar porque creo que la disrupcion tiene que caer en CDR3
    if disruption:
        productivity = 'potentially productive'
    else:
        productivity = 'productive'
        
    # get the sequences that we have to replace or insert and remove them in the consensus sequence
    clean_seq = gather_variants(list_insertions, fasta)
    clean_seq = gather_variants(list_deletions, clean_seq)
    # search for the motives and get CDR3
    ## modified sequence without indels
    CDR3, start, end, prod = get_CDR3(clean_seq, mincys)
    ## display also natural CDR3 without removing the variants
    CDR3_nature, start_nature, end_nature, prod_nature = get_CDR3(fasta, mincys)
    ## decide CDR3
    real_CDR3, real_start, real_end, prd = decide_CDR3((CDR3, CDR3_nature),
                                                       (start, start_nature), (end, end_nature),
                                                       (prod, prod_nature))
    if not prd:
        productivity = 'not productive'
    
    return real_CDR3, real_start, real_end, productivity

    # IGHD
    #fasta_IGHD = os.path.join('', fasta.replace('-fb.fa', '') + '_IGHD.fasta')
    #out_folder = blast_CDR3
    
    #if IGHD == '':
    #    # this IGHD sequence includes from CDR3 C to W in motif WGXG
    #    IGHD = clean_seq[real_start : real_end]
    
    #else:
    #    IGHD = IGHD
    
    #pipeline.CDR32fasta([IGHD], out_folder, fasta_IGHD)
    #blast_file = prepare_IGHD_blast(fasta_IGHD, out_folder)
    #print blast_file
    #dblast = blast_annotation(blast_file)

    # path for IMGT D fasta
    #emboss_list, EW_file = prepare_IGHD_emboss(fasta_IGHD, out_folder, '')
    
    #demboss = emboss_annotation(EW_file)

    #return real_CDR3, IGHD
    

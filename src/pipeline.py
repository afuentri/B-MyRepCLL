#! /usr/env/python                          #
# Afuentri                                  #
# Unidad de Genomica y Diagnostico Genetico #
# INCLIVA                                   #
# 2018-05-21                                #
#############################################

# Import modules
## python modules
import time
import sys
import argparse
import os
import logging as log
import string
import gzip
import re
import matplotlib
import glob
import glob2
import difflib
import subprocess
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Environment variables
repo = os.environ['scripts_repo']
## put references inside the repo
global IGH_references
IGH_references = os.environ['IGH_references']

## inner modules
global path_IGH_scripts
path_IGH_scripts = os.path.join(repo, 'pIgH/modules')
path_test = os.path.join(repo, 'pIgH/test')

# imports
import fastq_merge
import consensus2CDR3
from log_manager import log_setup
import arguments_parser_IGH
from dir_management import *

## Functions
def do_overwrite(folder, extension):

    """"""
    if owrt:
        overwrite = True
    else:
        overwrite = True
        files = os.listdir(folder)

        for f in files:
            if f.endswith(extension):
                overwrite = False
                break

    return overwrite


def execute(CMD):

    """"""
    subprocess.call([CMD], shell=True, executable='/bin/bash')


def create_CMD(name, folder):

    """"""
    name_cmd = str(cdate) + '-CMD_' + name + '.cmd'
    path_cmd = os.path.join(folder, name_cmd)
    name_log = str(cdate) + '-LOG_' + name + '.log'
    path_log = os.path.join(folder, name_log)

    return path_cmd, path_log


def process_communicate(CMD):

    """Return shell standard output"""
    proce = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE)
    
    return proce.communicate()[0].decode("utf-8").split(' ')


def removeEmpty(l):

    """Remove empty elements from list"""
    
    if '' in l:
        l.remove('')

    return l 


def write_header(fhand):

    """"""
    fhand.write('sample_name,Vregion,reads_mapped,reads_mapped_leader,reads_mapped_fr1,'
                'reads_mapped_fr2,reads_mapped_fr3,region_length,homology-IGHV_noFR3,mutational_status_noFR3,'
                                'joined,Jregion,J_assigned,J_coincidence,'
                                'IGHV-J,consensus_length,CDR3,CDR3_length,IGHD,IGHD_emboss,'
                                'insertions,deletions,ORF disruption,'
                'majorproductive_seq,nreads_majorproductive_seq,major_CDR3,major_IGHD,major_prod\n')


def dhom_write(s, v, k, d_hom):

    """Write hom output dict fields to tabular"""
    l = ('{},{},{},{},{},{},{},{},{},{},'
         '{},{},{},{},{},{},{},{},{},'
         '{},{},{},{},{},{},{},{},{}').format(s, v, d_hom[k]['nreads'], d_hom[k]['nreads_leader'],
                            d_hom[k]['nreads_fr1'], d_hom[k]['nreads_fr2'], d_hom[k]['nreads_fr3'],
                            d_hom[k]['length'], d_hom[k]['homology_noFR3'],
                            d_hom[k]['mutational_status_noFR3'],';'.join(d_hom[k]['joined']),
                            ';'.join(d_hom[k]['Jregion']), d_hom[k]['J_assigned'],
                            d_hom[k]['J_coincidence'],
                            d_hom[k]['IGHV-J'], str(d_hom[k]['consensus_length']),
                            d_hom[k]['CDR3'],
                            d_hom[k]['CDR3_length'], d_hom[k]['IGHD'],d_hom[k]['IGHD_emboss'],
                            d_hom[k]['insertions'],
                            d_hom[k]['deletions'], d_hom[k]['disruption'], d_hom[k]['new_seq'], d_hom[k]['nreads_new_seq'],
                            d_hom[k]['new_CDR3'], d_hom[k]['new_IGHD'], d_hom[k]['prod'])

    return l


def parse_vcf_IGHD(fpath, read, position=290):

    """"""
    l = read_file_simply_gzip(fpath)
    (list_insertions, list_deletions, disruption, 
     IGHD, pos_IGHD, pos_change, ins_IGHD) = [[], [], 0, '', '', 0, 0]

    for e in l:
        if not e.startswith("#"):
            nope, pos, nothing, ref, alt = e.split('\t')[:5]
            
            # look for indels
            if alt != "." and not ',' in alt:
                if len(alt) != len(ref):
                    if len(alt) > len(ref):
                        
                        # if the indel is an insertion bigger than 6 nucleotides and it is
                        # included in CDR3 region, define it as IGHD
                        if ((len(alt) - len(ref)) > 6) and alt in read and int(pos) > position:
                            IGHD = alt
                            pos_IGHD = pos
                            ins_IGHD = len(alt) - len(ref)
                        else:
                        
                            if ref in alt:
                                insertion = alt.replace(ref, '')

                                if len(insertion) % 3 != 0:
                                    disruption += 1
                            else:
                                if (len(alt) - len(ref)) % 3 != 0:
                                    disruption += 1

                            list_insertions.append(ref + '/' + alt + ':' + str(int(pos) + ins_IGHD))


                    elif len(alt) < len(ref):
                        if alt in ref:
                            deletion = ref.replace(alt, '')
                            if len(deletion) % 3 != 0:
                                disruption += 1
                        else:
                            if (len(ref) - len(alt)) % 3 !=0:
                                disruption += 1
                        list_deletions.append(ref + '/' + alt + ':' + str(int(pos) + ins_IGHD))

    if list_insertions == []:
        list_insertions = ['']
    if list_deletions == []:
        list_deletions = ['']
        
    return IGHD, pos_IGHD, list_insertions, list_deletions, disruption


def seqtk_quality_trimming(dic, out, merged_folder, CMD):

    """"""
    fcmd = open_file(CMD, mode='write')
    for fastq in dic:

        merged_name = os.path.join(merged_folder, os.path.basename(fastq))
        trimmed_name = os.path.join(out, dic[fastq]['trimmed_name'])
        CMD_trimming = 'seqtk trimfq {} | gzip > {}'.format(merged_name, trimmed_name)
        fcmd.write(CMD_trimming + '\n')

    fcmd.close()


def bbduk_trimming(dic, out, merged_folder, fastqs_folder, pairs, 
                   adapters_file, CMD, sides, mode='left'):

    """primer and quality trimming"""
    fcmd = open_file(CMD, mode='write')
    pairs_list = read_file_simply(pairs)
    trim_left = []
    for pair in pairs_list:
        p = pair.split(';')
        fastq1 = os.path.join(fastqs_folder, os.path.basename(p[0]))
        sample_name = os.path.basename(fastq1).split('_')[0]
        if mode == 'right':
            adapters_file = adapters_file.replace('left', 'right')
            
        merged_name = os.path.join(merged_folder, os.path.basename(fastq1))
        merged2_name = os.path.join(out, dic[fastq1]['trimmed_name']).replace('-trimmed', 
                                                                              '-trimmed5')
        if mode == 'left':
            if sides:
                trimmed_name = merged2_name
            else:
                trimmed_name = os.path.join(out, dic[fastq1]['trimmed_name'])
                
        elif mode == 'right':
            trimmed_name = os.path.join(out, dic[fastq1]['trimmed_name'])
            merged_name = os.path.join(merged_folder, os.path.basename(merged2_name))
            
        stats_file = trimmed_name.replace('.fastq.gz','-stats.txt')
        refstats_file = trimmed_name.replace('.fastq.gz','-refstats.txt')

        if len(p) > 1:
            # paired end mode
            fastq2 = os.path.join(fastqs_folder, os.path.basename(p[1]))
            merged_mate = os.path.join(merged_folder, os.path.basename(fastq2))
            merged2_mate = os.path.join(out, dic[fastq2]['trimmed_name']).replace('-trimmed', 
                                                                                  '-trimmed5')
    
            if mode == 'left':
                if sides:
                    trimmed_mate = merged2_mate
                else:
                    trimmed_mate = os.path.join(out, dic[fastq2]['trimmed_name'])
                    
                trim_left.append(trimmed_name + ';' + trimmed_mate)

            elif mode == 'right':

                trimmed_mate = os.path.join(out, dic[fastq2]['trimmed_name'])
                merged_mate =os.path.join(merged_folder, os.path.basename(merged2_mate))

            if mode == 'left':

                CMD_trimming = ('bbduk in1={} in2={} out1={} out2={} ref={} ktrim=l '
                                'k=21 hdist=1 mink=18 restrictleft=30 qtrim=rl'
                                ' trimq=15 minlength=50 stats={}'
                                ' refstats={}').format(merged_name, merged_mate, trimmed_name, 
                                                       trimmed_mate, adapters_file, stats_file, refstats_file)
            elif mode == 'right':

                CMD_trimming = ('bbduk in1={} in2={} out1={} out2={} ref={} ktrim=r'
                                ' k=21 hdist=1 mink=18 restrictright=30 qtrim=rl '
                                'trimq=15 minlength=50 stats={}'
                                ' refstats={}').format(merged_name, merged_mate, trimmed_name, 
                                                       trimmed_mate, adapters_file, stats_file, refstats_file)

        else:
            # single end mode
            if mode == 'left':
                CMD_trimming = ('bbduk in1={} out1={} ref={} ktrim=l k=21 hdist=1'
                                ' mink=18 restrictleft=30 qtrim=rl trimq=15 minlength=50'
                                ' stats={} refstats={}').format(merged_name, trimmed_name, 
                                                                adapters_file, stats_file, refstats_file)

            elif mode == 'right':
                CMD_trimming = ('bbduk in1={} out1={} ref={} ktrim=r k=21 hdist=1'
                                'mink=18 restrictright=30 qtrim=rl trimq=15 minlength=50'
                                'stats={} refstats={}').format(merged_name, 
                                                                            trimmed_name, 
                                                                            adapters_file, 
                                                                            stats_file, refstats_file)
            trim_left.append(trimmed_name)


        fcmd.write(CMD_trimming + '\n')

    fcmd.close()
    return trim_left


def flash(out_folder, trimmed_folder, pairs, CMD, new_pairs):

    """"""
    plist = read_file_simply(pairs)
    boole = True
    # check if they are paired end
    if len(plist[0].split(';')) > 1:
        fcmd = open_file(CMD, mode='write')
        fpairs = open_file(new_pairs, mode='write')
        for pair in plist:

            p = pair.split(';')
            fastq1_name = p[0]
            fastq2_name = p[1]
            fastq1_path = os.path.join(trimmed_folder, fastq1_name)
            fastq2_path = os.path.join(trimmed_folder, fastq2_name)

            CMD_flash = 'flash {} {} -t1 -M150 -z -o {} -d {}'.format(fastq1_path, fastq2_path,
                                                                      (fastq1_name.split('_')[0] + '_'), 
                                                                      out_folder)
            fcmd.write(CMD_flash + '\n')
            fpairs.write('%s\n' %(fastq1_name.split('_')[0] + '_.extendedFrags.fastq.gz'))
        fpairs.close()
        fcmd.close()

    else:
        # no need to merge because they are already single reads
        boole = False
        new_pairs = pairs
        fcmd = open_file(CMD, mode='write')
        fcmd.write('echo "Single end mode, no FLASH R1 and R2 merging"\n')
        fcmd.close()

    return new_pairs, boole


def bwamem_alignmentV(out, trimmed_folder, CMD, pairs, ref):

    """"""
    fcmd = open_file(CMD, mode='write')
    pairs_list = read_file_simply(pairs)
    for pair in pairs_list:
        p = pair.split(';')
        fastq1 = os.path.basename(p[0])
        sample_name = fastq1.split('_')[0]


        fastq1_path = os.path.join(trimmed_folder, fastq1)
        sortedbam_name = sample_name + '-sorted.bam'
        bam_path = os.path.join(out, sortedbam_name)
        tempbam_path = bam_path.replace('-sorted.bam', '-tempsorted.bam')
        sam = tempbam_path.replace('.bam', '.sam')
        if len(p) > 1:
            # paired end mode
            fastq2 = os.path.basename(p[1])
            fastq2_path = os.path.join(trimmed_folder, fastq2)
            #CMD_alignment = ("bwa mem -t4 -M -R "
            #                 "\"@RG\\tID:{0}\\tPL:ILLUMINA\\tSM:{0}\\tDS:ref={1}\\tCN:UGDG\\tDT:{5}\\tPU:{0}\" "
            #                 "{1} {2} {3} | samtools view -b - | samtools sort -o {4} -").format(sample_name, ref,
            #                                                                                      fastq1_path,
            #                                                                                      fastq2_path,
            #                                                                                      bam_path, cdate)
            
            ## filter FR3 reads
            CMD_alignment = ("bwa mem -t4 -M -R "
                             "\"@RG\\tID:{0}\\tPL:ILLUMINA\\tSM:{0}\\tDS:ref={1}\\tCN:UGDG\\tDT:{5}\\tPU:{0}\" "
                             "{1} {2} {3} | samtools view -b -o {7}; samtools view -H {7} > {6}; samtools view {7} | awk ' {{ if ($4 >= 200) print }} ' >> {6};"
                             "samtools view {7} | awk  ' {{ if ($4 < 200 && $6 !~ \"S\") print }} ' >> {6}; samtools view -b {6} | samtools sort -o {4} -; rm {6}; rm {7}").format(sample_name, ref,
                                                                                                 fastq1_path,
                                                                                                 fastq2_path,                                                                                                                                                                                                  bam_path, cdate,
                                                                                                 sam, tempbam_path)

        else:
            # single end mode
            #CMD_alignment = ("bwa mem -t4 -M -R "
            #                 "\"@RG\\tID:{0}\\tPL:ILLUMINA\\tSM:{0}\\tDS:ref={1}\\tCN:UGDG\\tDT:{4}\\tPU:{0}\" "
            #                 "{1} {2} |samtools view -b - | samtools sort -o {3} -").format(sample_name,
            #                                                                                ref, fastq1_path,
            #                                                                                bam_path, cdate)
            ## filter FR3 reads
            CMD_alignment = ("bwa mem -t4 -M -R "
                             "\"@RG\\tID:{0}\\tPL:ILLUMINA\\tSM:{0}\\tDS:ref={1}\\tCN:UGDG\\tDT:{4}\\tPU:{0}\" "
                             "{1} {2} | samtools view -b -o {6}; samtools view -H {6} > {5}; samtools view {6} | awk ' {{ if ($4 >= 200) print }} ' >> {5};"
                             "samtools view {6} | awk  ' {{ if ($4 < 200 && $6 !~ \"S\") print }} ' >> {5}; samtools view -b {5} | samtools sort -o {3} -; rm {5}; rm {6}").format(sample_name,
                                                                                             ref, fastq1_path, 
                                                                                             bam_path, cdate,
                                                                                             sam, tempbam_path)

        fcmd.write(CMD_alignment + '\n')

    fcmd.close()

def bwamem_alignmentJ(out, trimmed_folder, CMD, pairs, ref):

    """"""
    fcmd = open_file(CMD, mode='write')
    pairs_list = read_file_simply(pairs)
    for pair in pairs_list:
        p = pair.split(';')
        fastq1 = os.path.basename(p[0])
        sample_name = fastq1.split('_')[0]


        fastq1_path = os.path.join(trimmed_folder, fastq1)
        sortedbam_name = sample_name + '-sorted.bam'
        bam_path = os.path.join(out, sortedbam_name)

        if len(p) > 1:
            # paired end mode                                                                                                                                                                                  
            fastq2 = os.path.basename(p[1])
            fastq2_path = os.path.join(trimmed_folder, fastq2)
            CMD_alignment = ("bwa mem -t4 -M -R "
                             "\"@RG\\tID:{0}\\tPL:ILLUMINA\\tSM:{0}\\tDS:ref={1}\\tCN:UGDG\\tDT:{5}\\tPU:{0}\" "
                             "{1} {2} {3} | samtools view -b - | samtools sort -o {4} -").format(sample_name, ref,
                                                                                                 fastq1_path,
                                                                                                 fastq2_path,
                                                                                                 bam_path, cdate)

        else:
            # single end mode                                                                                                                                                                                  
            CMD_alignment = ("bwa mem -t4 -M -R "
                             "\"@RG\\tID:{0}\\tPL:ILLUMINA\\tSM:{0}\\tDS:ref={1}\\tCN:UGDG\\tDT:{4}\\tPU:{0}\" "
                             "{1} {2} | samtools view -b - | samtools sort -o {3} -").format(sample_name,
                                                                                             ref, fastq1_path,
                                                                                             bam_path, cdate)

        fcmd.write(CMD_alignment + '\n')

    fcmd.close()

    

def bams_list(out, pairs):

    """"""
    bams = []
    pairs_list = read_file_simply(pairs)
    for pair in pairs_list:
        p = pair.split(';')
        sample_name = p[0].split('_')[0]
        bam_name = sample_name + '-sorted.bam'
        bam_path = os.path.join(out, bam_name)
        bams.append(bam_path)

    return bams

def stat_list(bams, out):

    """"""
    stats = []
    for bam in bams:
        stats_name = os.path.basename(bam).replace('-sorted.bam', '-stats.txt')
        stats_path = os.path.join(out, stats_name)
        stats.append(stats_path)

    return stats


def bam_indexing(bams, CMD):

    """"""
    fcmd = open_file(CMD, mode='write')
    for bam in bams:
        CMD_indexing = 'samtools index {}'.format(bam)
        fcmd.write(CMD_indexing + '\n')
    fcmd.close()


def bam_stats(CMD, bams, out):

    """"""
    fcmd = open_file(CMD, mode='write')
    for bam in bams:
        stats_name = os.path.basename(bam).replace('-sorted.bam', '-stats.txt')
        stats_path = os.path.join(out, stats_name)
        CMD_stats = 'samtools idxstats {} > {}'.format(bam, stats_path)
        fcmd.write(CMD_stats + '\n')
    fcmd.close()


def list_vcfs(bams, folder):

    """"""
    vcf_list = []
    for bam in bams:
        #vcf = os.path.basename(bam).replace('-sorted.bam', '.vcf.gz')
        #vcf_path = os.path.join(folder, vcf)
        #vcf_list.append(vcf_path)
        fbvcf = os.path.basename(bam).replace('-sorted.bam', '-fb.vcf.gz')
        fbvcf_path = os.path.join(folder, fbvcf)
        vcf_list.append(fbvcf_path)

    return vcf_list


def parallel(CMD, p, LOG): ## make independent function with execute

    """"""
    parallel = ('parallel --compress --joblog {} '
                '-j {} --tmpdir {} :::: {}').format(LOG, p, TMPDIR, CMD)
    
    execute(parallel)

    
def parallel_resume(CMD, p, LOG):

    """"""
    parallel = ('parallel --compress --resume-failed --joblog {} '
                '-j {} --tmpdir {} :::: {}').format(LOG, p, TMPDIR, CMD)
    
    execute(parallel)


def list2file(l, fi, out):

    """"""
    fhand = os.path.join(out, fi)
    f = open_file(fhand, mode='write')
    f.write('\n'.join(l) + '\n')
    f.close()

    return fhand


def vcf_parsing(folder_results, vcf_folder, vcf_list, name):

    """"""
    # list coming from my file
    # remove not VCF files in VCF folder
    filename = name + '.csv'
    vcf_fof = list2file(vcf_list, 'vcf.fof', vcf_folder)
    out_file = os.path.join(folder_results, filename)
    script = os.path.join(path_IGH_scripts, 'samt-freeb_parser-onlyvariants.py')
    CMD = 'python {} -V {} -v > {}'.format(script, vcf_fof, out_file)
    execute(CMD)


def bcftools_consensus(fcmd, fasta_path, vcf, refV, current_ref, bam_path, sample):

    """"""
    if clonal:
        read_name = sample.split('_')[1]
        CMD_consensus = ('echo \">{3}\\n\" > {2}; samtools view {0}'
                         ' | grep {1} | cut -f10 >> {2}').format(bam_path, 
                                                                 read_name, 
                                                                 fasta_path, 
                                                                 current_ref)
    else:
        CMD_consensus = ('samtools faidx {} "{}" | bcftools'
                         ' consensus {} -o {}').format(refV, current_ref, 
                                                       vcf, fasta_path)
    fcmd.write(CMD_consensus + '\n')


def bcftools_consensus_complete(fcmd, bam, fasta_path_freebayes, path_current_ref, out):

    """"""
    vcf_temp_fb = os.path.basename(bam).replace('-sorted.bam', '-fb.vcf.gz')
    vcf_temp_fb_path = os.path.join(out, vcf_temp_fb)

    CMD_consensus_freebayes = ('freebayes -f {0} -m 0 --ploidy 1 --min-coverage {6}'
                               ' --genotype-qualities --strict-vcf -F {4} -C {5} {1} '
                               '| bgzip -c > {2}; tabix -f -p vcf {2} ; bcftools '
                               'consensus -f {0} {2} -o {3}').format(path_current_ref, 
                                                                     bam, vcf_temp_fb_path, 
                                                                     fasta_path_freebayes, 
                                                                     frec, minalt, mincoverage)
    fcmd.write(CMD_consensus_freebayes + '\n')

    return vcf_temp_fb_path


def prepare_trimming(out_folder, fastq_dictionary, merged_folder, 
                     fastqs_folder, ext, pairs, adapters, sides, 
                     prog='bbduk', mode='left'):

    """"""
    if len(ext) > 1:
        overwrite = do_overwrite(out_folder, '.fastq.gz')
        log.debug('The extensions found for FASTQ files are %s', list(ext))
    else:
        overwrite = do_overwrite(out_folder, list(ext)[0])
        log.debug('The extensions found for FASTQ files are %s', list(ext))
    if overwrite:
        log.debug('Starting FASTQ trimming. Out folder %s', out_folder)
        path_CMD_trimming, path_LOG_trimming = create_CMD('trimming', out_folder)
        if prog == 'seqtk':
            seqtk_quality_trimming(fastq_dictionary, out_folder, merged_folder, 
                                   path_CMD_trimming)
        elif prog == 'bbduk':
            if mode == 'left' and sides:
                path_CMD_trimming = path_CMD_trimming.replace('trimming', 'trimming5')
                path_LOG_trimming = path_LOG_trimming.replace('trimming', 'trimming5')
            
            list_trimmed = bbduk_trimming(fastq_dictionary, out_folder, merged_folder, 
                                          fastqs_folder, pairs, adapters, path_CMD_trimming, sides, mode)
        else:
            log.warning('No program named %s', prog)

        parallel(path_CMD_trimming, proc, path_LOG_trimming)
    else:
        log.warning('There are FASTQ files inside %s. Use mode overwrite if you '
                    'want to create FASTQ files in this folder', out_folder)


def prepare_flash(out, trimmed_folder, pairs):

    """"""
    new_pairs = os.path.join(out, 'pairs_trimming_flash.txt')
    overwrite = do_overwrite(out, '.fastq.gz')
    if overwrite:

        log.debug('Starting FASTQ read merging. Out folder %s', out)
        path_CMD_flash, path_LOG_flash = create_CMD('flash', out)
        new_pairs, boole = flash(out, trimmed_folder, pairs, path_CMD_flash, 
                                 new_pairs)
        parallel(path_CMD_flash, proc, path_LOG_flash)
    else:

        log.warning('There are FASTQ files inside %s. Use mode overwrite '
                    'if you want to create FASTQ files in this folder', out)
        fpairs = open_file(new_pairs, mode='write')
        p = glob.glob(str(out) + '/*.extendedFrags.fastq.gz')
        # if flash was not performed
        if len(p) == 0:
            p = glob.glob(str(trimmed_folder) + '/*.fastq.gz')
        for name in p:
            fpairs.write(os.path.basename(name) + '\n')

        fpairs.close()

    return new_pairs, boole


def split_fastqs(inf, out, pairs):

    """"""
    new_pairs_list = []
    pairs_list = read_file_simply(pairs)
    for l in pairs_list:
        name1 = l.split(';')[0]
        fastq1 = os.path.join(inf, name1)
        sample = name1.split('_')[0]
        rest = name1.split('_')[1:]
        if name1.endswith('.gz'):
            f1 = read_file_simply_gzip(fastq1)

        else:
            f1 = read_file_simply(fastq1)

        count = 0
        for s in range(0, len(f1), 4):
            count += 1
            f1_new = sample + '-' + str(count) + '_' + \
                     '_'.join(rest).replace('.gz','')
            f1_new_path = os.path.join(out, f1_new)
            fhand1 = open_file(f1_new_path, mode = 'write')
            read = f1[s:(s + 4)]
            fhand1.write('\n'.join(read))
            fhand1.close()

        if len(l.split(';')) > 1:

            name2 = l.split(';')[1]
            fastq2 = os.path.join(inf, name2)
            sample2 = name2.split('_')[0]
            rest2 = name2.split('_')[1:]
            if name2.endswith('.gz'):
                f2 = read_file_simply_gzip(fastq2)
            else:
                f2 = read_file_simply(fastq2)

            count2 = 0
            for m in range(0, len(f2), 4):
                count2 += 1
                f2_new = sample2 + '-' + str(count2) + '_' + \
                         '_'.join(rest2).replace('.gz','')
                f2_new_path = os.path.join(out, f2_new)
                fhand2 = open_file(f2_new_path, mode = 'write')
                read2 = f2[m:(m + 4)]
                pair = f1_new + ';' + f2_new
                new_pairs_list.append(pair)
                fhand2.write('\n'.join(read2))
                fhand2.close()
        else:
            pair = f1_new
            new_pairs_list.append(pair)

    return new_pairs_list


def clonality_trimmed(trimmed_folder, out, pairs):

    """"""
    overwrite = do_overwrite(out, '.fastq.gz')
    pairs_trimming = os.path.join(out, 'pairs_trimming.txt')
    fpairs_trimming = open_file(pairs_trimming, mode='write')
    if overwrite:
        log.debug('Trimmed FASTQ preparation for clonality analysis. Out folder %s', out)
        new_list_pairs = split_fastqs(trimmed_folder, out, pairs)
        fpairs_trimming.write('\n'.join(new_list_pairs))
        fpairs_trimming.close()
    else:
        log.warning('There are FASTQ files inside %s. Use mode overwrite if '
                    'you want to create FASTQ files in this folder', out)
        p = glob.glob(str(out) + '/*')
        for name in p:
            pair_exists = False
            if os.path.basename(name).endswith('.fq') or os.path.basename(name).endswith('.fq.gz') or\
               os.path.basename(name).endswith('.fastq') or os.path.basename(name).endswith('.fastq.gz'):
                if ('_R1' in name) or ('_1' in name):
                    for name2 in p:
                        if name.split('_')[0] == name2.split('_')[0]:
                            pair_exists = True
                            fpairs_trimming.write(os.path.basename(name) + ';' + os.path.basename(name2))
                    if not pair_exists:
                        fpairs_trimming.write(os.path.basename(name))


    return pairs_trimming


def prepare_alignment(out_folder, trimmed_folder, pairs, ref):

    """"""
    # BWA mem can have inner parallelization
    overwrite = do_overwrite(out_folder, '.bam')
    if overwrite:
        log.debug('Starting FASTQ alignment. Out folder %s', out_folder)
        path_CMD_alignment, path_LOG_alignment = create_CMD('alignment', out_folder)
        bam_list = bams_list(out_folder, pairs)
        if 'bamsV' in out_folder:
            bwamem_alignmentV(out_folder, trimmed_folder, path_CMD_alignment, pairs, ref)
        elif 'bamsJ' in out_folder:
            bwamem_alignmentJ(out_folder, trimmed_folder, path_CMD_alignment, pairs, ref)
        p = int(int(proc)/4) # consider inner threading bwa mem
        parallel(path_CMD_alignment, p, path_LOG_alignment)
        path_CMD_bamindex, path_LOG_bamindex = create_CMD('bamindex', out_folder)
        bam_indexing(bam_list, path_CMD_bamindex)
        parallel(path_CMD_bamindex, proc, path_LOG_bamindex)

    else:
        bam_list = bams_list(out_folder, pairs)
        log.warning('There are BAM files inside %s. Use mode overwrite '
                    'if you want to create BAM files in this folder', out_folder)

    return bam_list


def prepare_stats(out_folder, fastq_dictionary, alignment_folder, bam_list):

    """"""
    # BAMS HAVE TO BE SORTED
    overwrite = do_overwrite(out_folder, '-stats.txt')
    if overwrite:
        log.debug('Starting mapping statistics from sorted BAM files. '
                  'Out folder %s', out_folder)
        path_CMD_stats, path_LOG_stats = create_CMD('stats', out_folder)
        stats_list = stat_list(bam_list, out_folder)
        bam_stats(path_CMD_stats, bam_list, out_folder)
        parallel(path_CMD_stats, proc, path_LOG_stats)
    else:
        stats_list = stat_list(bam_list, out_folder)
        log.warning('There are stats files inside %s. Use mode overwrite if '
                    'you want to create stats files in this folder', out_folder)

    return stats_list


def prepare_subtype(out_folder, list_bamp):

    """"""
    overwrite = do_overwrite(out_folder, '.txt')
    if overwrite:
        log.debug('Starting subtype resume for bam parsing. Out folder %s', out_folder)
        path_CMD_resume, path_LOG_resume = create_CMD('resume', out_folder)
        bam_parsing_resume(path_CMD_resume, list_bamp, out_folder)
        parallel(path_CMD_resume, proc, path_LOG_resume)
    else:
        log.warning('There are output files inside %s. Use mode overwrite if you want'
                    ' to create files in this folder', out_folder)

    subtype_list = []
    subtype_list = glob.glob(str(out_folder) + '/subtype*')

    return subtype_list


def prepare_joined_references(out_folder, ref1, ref2):

    """"""
    overwrite = do_overwrite(out_folder, 'IMGT-IGHV-J.fa')
    if overwrite:
        log.debug('Starting generation of the joined references of '
                  'IGHV and IGHJ. Out folder %s', out_folder)
        out_name = os.path.join(out_folder, 'IMGT-IGHV-J.fa')
        f1 = read_file_simply(ref1)
        f2 = read_file_simply(ref2)
        f_out = open_file(out_name, mode = 'write')
        for line in f1:
            if line.startswith('>'):
                shortV = line.replace('>','').replace('/','-')

            else:
                seq1 = line.strip()
                for line2 in f2:

                    if line2.startswith('>'):
                        shortJ = line2.replace('>','')
                        f_out.write('>' + shortV + '_' + shortJ + '\n')

                    else:
                        seq2 = line2.strip()
                        f_out.write(seq1 + seq2 + '\n')
        f_out.close()
        index_ref = 'bwa index {}'.format(out_name)
        execute(index_ref)

    else:
        out_name = os.path.join(out_folder, 'IMGT-IGHV-J.fa')
        log.warning('There is already a FASTA file called IMGT-IGHV-J.fa in %s. '
                    'Use mode overwrite if you want to overwrite the existing file.', out_folder)


    return out_name


def probable_regions(stats, out, name):

    """"""
    fname = os.path.join(out, name)
    script = os.path.join(path_IGH_scripts, 'probable_regions.py')
    CMD = 'python {} -f {} -v -o {}'.format(script, stats, fname)
    log.info('Computing probable regions for %s', stats)
    execute(CMD)

    return fname


def bam_parsing(fof, out):

    """"""
    overwrite = do_overwrite(out, '.txt')
    if overwrite:
        log.debug('Starting bam parsing. Out folder is %s', out)
        script = os.path.join(path_IGH_scripts, 'IGHBamsParser_v2.py')
        if clonal:
            CMD = 'python {} -V {} -v -o {} --clonal'.format(script, fof, out)
            log.info('Computing bam parsing in clonality mode for %s', fof)
        else:
            CMD = 'python {} -V {} -v -o {}'.format(script, fof, out)
            log.info('Computing bam parsing for %s', fof)
        execute(CMD)

    else:
        log.warning('There are already bam parsing output files in %s. '
                    'Use mode overwrite if you want to overwrite the existing file.', out)

    # obtain a list with output files
    list_bampars = glob.glob(str(out) + '/info*')

    return list_bampars


def bam_parsing_resume(CMD, list_pars, out):

    """"""
    fcmd = open_file(CMD, mode='write')
    for l in list_pars:
        l_path = l
        outname = os.path.basename(l).replace('info_bams_', 'subtypes_resume_')
        outpath = os.path.join(out, outname)
        CMD_resume = ('cat {} | grep "region" | cut -f2,4 | sort | uniq -c | sort'
                      ' -rn > {}').format(l_path, outpath)
        fcmd.write(CMD_resume + '\n')

    fcmd.close()


def pie_plots(d, out):

    """"""
    overwrite = do_overwrite(out, '.png')
    if overwrite:
        log.debug('Starting pie plot generation. Out folder is %s', out)
        for sample in d:
            if d[sample] != {}:
                # get keys as labels in the order of numerically ordered values
                l = sorted(d[sample], key=lambda i: float(d[sample][i]))
                p = sorted(d[sample].values())
                # Define the basal percent which is the rest to 100%
                rest = float(100-sum(p))
                l.append('basal')
                p.append(rest)
                # make a pie plot for each sample with the represented V genes
                #labels = ','.join(l)

                ##  highlight the highest percent
                wid = 0
                explode = ((wid,) * (len(l) - 2)) + (0.05, 0)

                ## Save the plot in the given out folder
                name = sample + '_representedVgenes.png'
                path_name = os.path.join(out, name)
                plt.pie(p, explode=explode, labels=l, autopct='%1.1f%%', shadow=True)
                plt.savefig(path_name)
                plt.gcf().clear()
    else:
        log.warning('There are already pie plot files in %s. '
                    'Use mode overwrite if you want to overwrite the existing file.', out)

    # obtain a list with output files
    pies = glob.glob(str(out) + '/*_representedVgenes.png')

    return pies


def VJref_dictionary(ref):

    """Dictionary containing IGHV alleles"""
    # Open file and store it to list
    ref_list = read_file_simply(ref)

    # Looping through list and storing elements in dictionary
    dict = {}
    for element in ref_list:
        element = element.strip()
        if element.startswith('>'):
            # ref name will be key
            cabecera = element.replace('>','').replace('/','-')
        else:
            # concatenate if there is more than one sequence line per reference
            if cabecera not in dict:
                dict[cabecera] = element
            else:
                dict[cabecera] = str(dict[cabecera]) + element
            # in Vref_dictionary else statement was avoided

    return dict


def make_reference(path, cabecera, sequence, mode):

    """Write reference sequence for alignment to file"""
    f = open_file(path, mode)
    f.write('>' + cabecera + '\n' + sequence + '\n')
    f.close()

    return f


def calculate_basal(bam_listV, regions_listV):

    """
    If given a policlonal control, replace the default basal 
    percent by a new calculated basal percent
    """
    nb = []
    new_basal = basal
    policlonal = False
    log.info('Default basal percent of reads is %s', new_basal)
    for bam in bam_listV:
        if 'policlonal' in bam:
            policlonal = True
            sample_name = os.path.basename(bam).replace('-sorted.bam', '')
            log.info('policlonal sample detected: %s, calculating new basal filter...', 
                     bam)
            total_reads = 0
            d = {}
            percent_list = []

            for reg in regions_listV:
                if reg.split(',')[0] == sample_name:
                    s, ref, counts, length = reg.split(',')
                    Vgene = '-'.join(ref.split('_')[0].split('-')[:-1])
                    total_reads += int(counts)
                    if not Vgene in d:
                        d[Vgene] = int(counts)
                    else:
                        d[Vgene] = d[Vgene] + int(counts)

            for k in d:
                percent_list.append(float(d[k])/float(total_reads)*100)
            if percent_list != []:
                nb.append(max(percent_list))


    if policlonal:
        new_basal = max(nb)
        log.info('Applying new basal filter: %s percent', new_basal)


    return new_basal                                


def filter_Vregion_basal(regions_listV, bam_listV):

    """
    Filter the list with the represented regions of IGHV ignoring 
    a basal percent of reads from policlonal control
    """
    c = {}
    max_allele = {}

    if basal != 0:
        nbasal = calculate_basal(bam_listV, regions_listV)
    else:
        nbasal = basal

    for bam in bam_listV:
        sample_name = os.path.basename(bam).replace('-sorted.bam', '')

        if not sample_name in c:
            c[sample_name] = {}
            max_allele[sample_name] = {}
            c[sample_name]['total_counts'] = 0
        for reg in regions_listV:
            if reg.split(',')[0] == sample_name:
                s, ref, counts, length, leader, fr1, fr2, fr3 = reg.split(',')
                ref = ref.replace('D-','-D') #include D name in the allele field
                Vgene = '-'.join(ref.split('_')[0].split('-')[:-1])

                if not Vgene in c[sample_name]:
                    # the option is this or a dictionary per sample
                    max_allele[sample_name][Vgene] = ref + ';' + counts
                    c[sample_name][Vgene] = int(counts)
                else:
                    c[sample_name][Vgene] = c[sample_name][Vgene] + int(counts)
                    if int(max_allele[sample_name][Vgene].split(';')[1]) > int(counts):
                        max_allele[sample_name][Vgene] = max_allele[sample_name][Vgene]
                    else:
                        max_allele[sample_name][Vgene] = ref + ';' + counts

                c[sample_name]['total_counts'] += int(counts)

    new, new_gene, c_counts, c_counts_nofilter = ([], [], {}, {})
    
    for bam in bam_listV:
        sample_name = os.path.basename(bam).replace('-sorted.bam', '')
        if not sample_name in c_counts:
            c_counts[sample_name] = {}
        if not sample_name in c_counts_nofilter:
            c_counts_nofilter[sample_name] = {}
        for reg in regions_listV:
            if reg.split(',')[0] == sample_name:

                s, ref, counts, length, leader, fr1, fr2, fr3 = reg.split(',')
                ref = ref.replace('D-','-D') #include D name in the allele field
                Vgene = '-'.join(ref.split('_')[0].split('-')[:-1])
                if max_allele[sample_name][Vgene].split(';')[0] == ref:
                    counts = c[sample_name][Vgene]
                    new_gene.append('%s,%s,%s,%s,%s,%s,%s,%s' %(s, ref.replace('-D','D-'), counts, length,
                                                                leader, fr1, fr2, fr3))
                perc = (float(c[sample_name][Vgene])/float(c[sample_name]['total_counts']))*100
                if not Vgene in c_counts_nofilter[sample_name]:
                    c_counts_nofilter[sample_name][Vgene] = perc

                if perc >= nbasal:
                    new.append(reg)
                    if not Vgene in c_counts[sample_name]:
                        c_counts[sample_name][Vgene] = perc



    log.info('Dictionary total counts: %s', c)
    return new, new_gene, c_counts, c_counts_nofilter, nbasal


def filter_Jregion_basal(regions_listJ, bam_listV):

    """Filter the list with the represented regions of IGHJ ignoring a basal
 percent of reads from policlonal control"""
    c = {}
    for bam in bam_listV:
        sample_name = os.path.basename(bam).replace('-sorted.bam', '')
        for reg in regions_listJ:

            if reg.split(',')[0] == sample_name:

                s, ref, counts, leader, fr1, fr2, fr3 = reg.split(',')
                if not sample_name in c:
                    c[sample_name] = int(counts)
                else:
                    c[sample_name] = c[sample_name] + int(counts)

    new = []

    for bam in bam_listV:
        sample_name = os.path.basename(bam).replace('-sorted.bam', '')

        if sample_name in c:
            if c[sample_name] > 0:
                for reg in regions_listJ:
                    if reg.split(',')[0] == sample_name:
                        s, ref, counts, leader, fr1, fr2, fr3 = reg.split(',')
                        perc = (float(counts)/float(c[sample_name]))*100

                        if perc >= basal:
                            new.append(reg)


    log.info('Dictionary total counts J: %s', c)
    return new


def d_output_creation():

    d = {'length' : '', 'nreads' : '', 'nreads_leader' : '',
         'nreads_fr1' : '', 'nreads_fr2' : '', 'nreads_fr3' : '',
         'homology_freebayes' : '',
         'mutational_status_freebayes' : '',
         'joined' : '', 'Jregion' : '', 'J_assigned' : '',
         'J_coincidence' : '', 'IGHV-J': '',
         'consensus_length': '', 'CDR3_length' : '',
         'CDR3' : '', 'IGHD' : '', 'IGHD_emboss' : '',
         'insertions':'', 'deletions':'', 'disruption':'',
         'ori_junction':'', 'perc_ori_junction':'', 'new_seq':'',
         'nreads_new_seq':'', 'new_CDR3':'', 'new_IGHD':'',
         'prod':''}

    return d


def annotate_hom(bam_listV, regions_listJ, regions_listV, d_hom):

    """"""
    for bam in bam_listV:
        sample_name = os.path.basename(bam).replace('-sorted.bam', '')
        if not clonal:
            J_region = J_region_annotation(regions_listJ, sample_name)
        for reg in regions_listV:
            if reg.split(',')[0].split('_')[0] == sample_name: ##clonal
                s, ref, counts, length, leader, fr1, fr2, fr3 = reg.split(',')
                
                if clonal:
                    J_region = J_region_annotation(regions_listJ, s)
                # dictionary key comprised by sample name and Vregion
                sample_region = s.replace(':','-') + '_' + ref # clonal
                # create the subkeys for the previous key
                d_hom[sample_region] = d_output_creation()

                d_hom[sample_region]['length'] = length
                d_hom[sample_region]['nreads'] = counts
                d_hom[sample_region]['nreads_leader'] = leader
                d_hom[sample_region]['nreads_fr1'] = fr1
                d_hom[sample_region]['nreads_fr2'] = fr2
                d_hom[sample_region]['nreads_fr3'] = fr3
                d_hom[sample_region]['Jregion'] = J_region

    return d_hom


def mutationalStatus(consensus_complete_list, homology_folder, ref_dict):

    """"""
    overwrite = do_overwrite(homology_folder, '.water')
    water_list = []
    if overwrite:
        log.debug('Starting IGHV mutational status characterization. '
                  'Out folder is %s', homology_folder)
        path_CMD_homology, path_LOG_homology = create_CMD('homology',
                                                          homology_folder)
        fcmd_homology = open_file(path_CMD_homology, mode='write')

        for consensus_path_fb in consensus_complete_list:
            if not 'gene' in os.path.basename(consensus_path_fb):
                sample_name = os.path.basename(consensus_path_fb).split('_')[0]

                # local alignment of consensus sequences with the references
                water_out_fb = homology(ref_dict, consensus_path_fb,
                                        fcmd_homology, sample_name, homology_folder)

                ## add local alignment results files to list
                water_list.append(water_out_fb)

        fcmd_homology.close()
        log.info('Executing %s; calculating homology for each sample'
                 ' and probable allele', path_CMD_homology)
        parallel(path_CMD_homology, proc, path_LOG_homology)

    else:
        log.warning('There are already homology output files in %s. '
                    'Use mode overwrite if you want to overwrite the '
                    'existing file.', homology_folder)

        water_list = glob2.glob(homology_folder + '/**/*.water' )
    
                

    return water_list


def J_region_annotation(regions_listJ, sample):

    """"""
    J = {}
    n = []
    J_region = []
    for m in range(len(regions_listJ)):
        l = regions_listJ[m]

        if l.split(',')[0] == sample: # clonal
            name, jref, num, leader, fr1, fr2, fr3 = l.split(',')
            if int(num) > 0:
                J[num] = jref
                n.append(int(num))

    for i in n:

        J_region.append(str(i) + ' ' + J[str(i)])

    if not J_region:
        J_region = ''
        
    return J_region


def bam2sam(bam_path, ref, sams_path):

    """"""
    CMD = 'samtools view -h {}  "{}" > {}'.format(str(bam_path), str(ref), str(sams_path))
    execute(cmd_regions)


def consensus_sequenceV(ref, ref_dict, out, bam_path, vcf_folder, refV, sample, fcmd):

    """"""
    # Reference for V region extraction
    ## ref is the original ref name found in BAM files
    ## reference allele nuc sequence
    short_ref = ref.replace('/','-')
    current_ref = ref_dict[short_ref]

    # Paths and commands
    ## consensus files
    fasta_consensus_name_freebayes = sample.replace(':','-') + '_' + short_ref + '-fb.fa'
    fasta_path_freebayes = os.path.join(out, fasta_consensus_name_freebayes)

    # Consensus sequence
    vcf_fb = os.path.basename(bam_path).replace('-sorted.bam', '-fb.vcf.gz')
    vcf_fb_path = os.path.join(vcf_folder, vcf_fb)
    bcftools_consensus(fcmd, fasta_path_freebayes, vcf_fb_path,  refV, ref, bam_path, sample)

    return fasta_path_freebayes


def consensus_sequence_complete(ref, ref2, out, folder_references, bam, sample_name, fcmd, out_vcfs):

    """"""
    # Reference extraction
    ref_path = os.path.join(folder_references, (ref + '.fa'))

    # Paths and commands
    ## consensus files
    fasta_consensus_name_freebayes = sample_name + '_' + ref2 + '-fb.fa'
    fasta_path_freebayes = os.path.join(out, fasta_consensus_name_freebayes)

    # Consensus sequence
    vcf_path = bcftools_consensus_complete(fcmd, bam, fasta_path_freebayes, ref_path, out_vcfs)
    
    return fasta_path_freebayes, vcf_path


def homology(ref_dict, consensus_path, fcmd_homology, sample_name, homology_folder):

    """"""
    # Create a directory for each sample inside homology
    dir_sample = os.path.join(homology_folder, sample_name)
    create_dir(dir_sample)
    # allele sequence extraction of the V dictionary
    ## clonality mode key ID is sample_read_Vgene
    if clonal:
        short_ref = os.path.basename(consensus_path).split('_')[2]
    else:
        short_ref = os.path.basename(consensus_path).split('_')[1]
        
    if 'gene' in short_ref:
        short_ref = short_ref.replace('gene', '01')
        
    current_ref = ref_dict[short_ref]
    reference_name = short_ref + '.fa'
    reference_path = os.path.join(homology_folder, reference_name)
    out_water = os.path.basename(consensus_path).replace('.fa','.water')
    out_path = os.path.join(dir_sample, out_water)
    if not os.path.isfile(reference_path):
        make_reference(reference_path, short_ref, current_ref, mode='write')

    cmd_emboss = ('emboss_water -asequence {} -bsequence {} -gapopen 10.0 '
                  '-gapextend 0.5 -outfile {} '
                  '-aformat pair').format(reference_path, consensus_path, out_path)
    fcmd_homology.write(cmd_emboss + '\n')

    return out_path



def homology_parsing(alignment_list, d, g, field, field2):

    """Parse EMBOSS water format 'pair' from a list of paths"""
    for path in alignment_list:
        path = path.strip()
        l = read_file_simply(path)
        name = os.path.basename(path)
        if name.endswith('-fb.water'):
            k = '_'.join(name.replace('-fb.water', '').split('_')[:2])
            if clonal:
                
                k = '{}_{}'.format(k.split('_')[0].replace('123456789',
                                                           '-'), k.split('_')[1])
                
        
     
        if k in d:
        
            for e in range(len(l)):
                if l[e].startswith('# Length'):
                    h = l[e].replace("#","") + ';' + l[e + 1].replace("#","")
                    hom = h.split('(')[1].replace(")","").replace("%","")
                    if float(hom) > 85:
                        a = '(valid)'
                        if float(hom) >= 98:
                            mut_status = 'UM'
                        else:
                            mut_status = 'MM'
                    else:
                        a = '(not valid)'
                        mut_status = ''


            h = h + a
            
            # Write homology value to dictionary
            d[k][field] = h
            d[k][field2] = mut_status
            if k in g:
                g[k][field] = h
                g[k][field2] = mut_status

    ## include empty homology for samples in dictionary without homology data
    for r in d.keys():
        if not field in d[r]:
            d[r][field] = ''
        if not field2 in d[r]:
            d[r][field2] = ''
    
    for r in g.keys():
        if not field in g[r]:
            g[r][field] = ''
        if not field2 in g[r]:
            g[r][field2] = ''

    return d, g


def joined_annotation(d_hom, g_hom, path_list):

    """"""
    for path in path_list:
        sample = os.path.basename(path).replace('subtypes_resume_','').replace('.txt','')

        for k in d_hom:
            if k.split('_')[0] == sample:
                
                reg = k.split('_')[1]

                CMD = 'grep "{}" {} | grep -v "None"'.format(reg, path)
                proce = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE)
                out = proce.communicate()[0].decode("utf-8").split('\n')
                
                l = ''

                n = removeEmpty(out)
                
                if len(d_hom[k]['Jregion']) > 0:
                    # check if not Jregion identified
                    if len(n) == 0:
                        
                        d_hom[k]['J_assigned'] = d_hom[k]['Jregion'][0].split(' ')[1]
                        d_hom[k]['J_coincidence'] = 'no joined info'
                        
                    else:
                        d = {}
                        for e in range(len(n)):

                            n[e] = n[e].replace('\r','').replace('\e','').strip()
                            m = n[e].split('\t')[1].split(';')[0]
                            co = int(n[e].split('\t')[0].split(' ')[0])
                            if not m in d:
                                d[m] = co
                            else:
                                d[m] += co

                        ## pick maximum JH representation for that allele
                        ma = max(d, key=d.get)
                    
                        min_joined_reads = float(d_hom[k]['nreads'])*0.8/100
                        if d[ma] >= min_joined_reads:
                                                                       
                            d_hom[k]['J_assigned'] = ma
                            # assign IGHJ and determine the coincidence with the
                            # joined subtype (overlapping reads between V and J)
                            if d_hom[k]['J_assigned'] == d_hom[k]['Jregion'][0].split(' ')[1]:
                                d_hom[k]['J_coincidence'] = 'yes'
                            else:
                                d_hom[k]['J_coincidence'] = 'no'

                            l = '{} {}'.format(str(d[ma]), ma)

                        else:

                            d_hom[k]['J_assigned'] = d_hom[k]['Jregion'][0].split(' ')[1]
                            d_hom[k]['J_coincidence'] = 'not passing filters'
                        
                else:
                    d_hom[k]['J_assigned'] = ''
                    d_hom[k]['J_coincidence'] = ''
                d_hom[k]['joined'] = [l]
                regV = k.split('_')[-1]
                d_hom[k]['IGHV-J'] = regV + '_' + d_hom[k]['J_assigned']
                if k in g_hom:
                    g_hom[k]['IGHV-J'] = d_hom[k]['IGHV-J']
                    g_hom[k]['joined'] = d_hom[k]['joined']
                    g_hom[k]['J_assigned'] = d_hom[k]['J_assigned']
                    g_hom[k]['J_coincidence'] = d_hom[k]['J_coincidence']
    
    return d_hom, g_hom


def fastq2fasta(fastq, fasta_path, sample_name):

    """"""
    fasta = os.path.join(fasta_path, (sample_name + '.fa'))
    ffasta = open_file(fasta, mode='write')
    m = read_file_simply_gzip(fastq)
    seq = m[1]
    ffasta.write('>seq1' + '\n' + seq + '\n')
    ffasta.close()

    return fasta


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


def local_alignmentJ(fastq, out, refJ, sample_name):

    """"""
    identity_dict = {}
    all_identities = {}
    scores = {}
    fasta = fastq2fasta(fastq, out, sample_name)
    resume = fasta.replace('.fa', '_resume.txt')
    out_resume = open_file(resume, mode='write')
    out_ref = os.path.join(out, 'ref.fa')
    out_path = os.path.join(out, 'out_emboss.txt')
    reference_path = refJ
    ref = read_file_simply(refJ)
    for read in range(len(ref)):
        refseq = open_file(out_ref, mode='write')
        if ref[read].startswith('>'):
            ref_name = ref[read].replace('>', '')
            refseq.write(ref[read] + '\n' + ref[read + 1] + '\n')
            refseq.close()

            CMD = ('emboss_water -asequence {} -bsequence {} -gapopen 10.0 '
                   '-gapextend 0.5 -outfile {} '
                   '-aformat pair').format(out_ref, fasta, out_path)
            execute(CMD)

            identity, whole, score = parse_emboss(out_path)
            all_identities[ref_name] = float(identity)
            scores[ref_name] = float(score)
            if float(score) > 90.0 and float(identity) > 60.0:
                identity_dict[ref_name] = float(identity)
            out_resume.write(ref_name + ':' + whole + ';' + score + '\n')

    if identity_dict == {}:
        max_value = max(all_identities.values())
        max_identity = [key for key in all_identities.keys() if all_identities[key] == max_value]
    else:
        max_value = max(identity_dict.values())
        max_identity = [key for key in identity_dict.keys() if identity_dict[key] == max_value]

    # subset scores dictionary to the IGHJ alleles with max value of homology
    max_scores = dict((k, scores[k]) for k in max_identity if k in scores)
    # choose the biggest score when identities are the same
    max_identity_J = max(max_scores, key=max_scores.get)
    out_resume.close()


def alternative_J_assign(d_hom, pairs, out, trimmed_folder, refJ):

    """Local alignment algorithm assignation of J region"""
    for s in d_hom:
    
        if d_hom[s]['J_assigned'] == '':
            
            Jallele = ''
            Jallele2 = ''
            pairs_list = read_file_simply(pairs)
            for pair in pairs_list:
                p = pair.split(';')
                fastq1 = os.path.basename(p[0])
                sample_name = fastq1.split('_')[0]
                if sample_name == s.split('_')[0]:
                    fastq1_path = os.path.join(trimmed_folder, fastq1)
                    m = read_file_simply_gzip(fastq1_path)
                    if len(m) == 4:
                        Jallele = local_alignmentJ(fastq1_path, out, refJ, s)

                        if len(pair.split(';')) > 1:
                            fastq2 = os.path.basename(p[1])
                            fastq2_path = os.path.join(trimmed_folder, fastq2)
                            Jallele2 = local_alignmentJ(fastq2_path, out, refJ, s)

            if not Jallele:
                Jallele = ''
            d_hom[s]['J_assigned'] = Jallele
            d_hom[s]['Jregion'] = [Jallele, Jallele2]
            d_hom[s]['IGHV-J'] = '{}_{}'.format(s.split('_')[1], d_hom[s]['J_assigned'])

    return d_hom


def CDR32fasta(list_reads, out, fasta_path):

    """"""
    print(list_reads)
    fasta = open_file(fasta_path, mode = 'write')
    for i in range(len(list_reads)):
        if list_reads[i] != '':
            fasta.write('>' + str(i + 1) + '\n' + list_reads[i] + '\n')


def info_motifs(fhand, vcf, d_hom, g_hom):

    """"""
    file_name = os.path.basename(fhand)
    vcf_name = file_name.replace('.fa','.vcf.gz').replace('-sm','')
    vcf_complete_path = os.path.join(vcf, vcf_name)
    file_data = file_name.replace('-fb.fa','')
    file_name2 = file_name.replace('.fa', '')

    if clonal:
        sample_name = '_'.join(file_data.split('_')[0:2])
        Vallele = file_data.split('_')[2]
        
        k = sample_name.replace('123456789', '-') + '_' + Vallele
        
        Jallele = d_hom[k]['J_assigned']
        genes = False
        # complete
    else:
        sample_name = file_data.split('_')[0]
        Vallele = file_data.split('_')[1]
        Jallele = file_data.split('_')[2]
        VJallele = file_data.split('_')[1:]
        if 'gene' in Vallele:
            genes = True
            for e in g_hom:
                if  e.split('_')[0] == sample_name:
                    if (Vallele.replace('gene','')) in e.replace('D-','-D'):
                        k = e
        else:
            genes = False
            k = '_'.join(file_data.split('_')[:2])

    return file_name, file_name2, Vallele, vcf_complete_path, k, genes


def function_motifs(k, f, file_name2, out, Vallele, vcf_complete_path):

    """"""
    # definition of empty lists
    read, IGHD = ('', '')
    flist = read_file_simply(f)
    info_file = file_name2 + '_infoCDR3.txt'
    info_path = os.path.join(out, info_file)
    finfo = open_file(info_path, mode='write')
    fasta_name = file_name2 + '.fasta'
    fastaD_name = file_name2 + '-IGHD.fasta'
    #fastaDaa_name = file_name2 + '-IGHD.fasta'
    fasta_path = os.path.join(out, fasta_name)
    fastaD_path = os.path.join(out, fastaD_name)
    
    # Get sequence and clean it
    # function
    s = []
    for m in fist[1:]:
        if not m.startswith('>'):
            s.append(m.strip())
        else:
            break

    # isolate sequence to extract CDR3
    sequence = ''.join(s)
    sequence = sequence.upper().strip()

    # get CDR3 sequence
    aa, start, end, prod = consensus2CDR3.cdr3_extraction(sequence, mincys=3)
    
    if aa:
        read = sequence[start:end]
        
    if clonal:
        # we do not have VCF
        list_ins, list_dels, disr = ([], [], '')
        if aa:
            IGHD = sequence[(start + 9) : (end - 18)]
            
    else:
        IGHD, pos_IGHD, list_ins, list_dels, disr = parse_vcf_IGHD(vcf_complete_path,
                                                                   sequence, 290)
    if IGHD == '':
        IGHD = calculateD(sequence, start, end)
        

    finfo.close()
    
    return (read, IGHD, aa, fasta_path, fastaD_path,
            list_ins, list_dels, disr)


def motifs(seq_list, d_hom, g_hom, out, vcf_path):

    """"""
    #seq_list is a fasta file
    fasta_CDR3_list = []
    fasta_IGHD_list = []
    
    junction_seqs = {}
    ## iterate over given sequences to extract info
    for f in seq_list:
        
        (file_name, file_name2, Vallele,
         vcf_complete_path, k, genes) = info_motifs(f, vcf_path,
                                                    d_hom, g_hom)
        # close info_motifs
        # check if key in d_hom
        if k in d_hom:
            
            (read, IGHD, aa, fasta_path,
             fastaD_path, list_ins, list_dels, disr) = function_motifs(k, f,
                                                                       file_name2, out,
                                                                       Vallele, 
                                                                       vcf_complete_path)
            junction_seqs[k] = read
            CDR32fasta([read], out, fasta_path)
            CDR32fasta([IGHD], out, fastaD_path)
            fasta_IGHD_list.append(fastaD_path)
            fasta_CDR3_list.append(fasta_path)

            if 'fb' in file_name:
                field, field2, field3, field4, field5 = (['CDR3', 'CDR3_length', 
                                                         'insertions', 'deletions',
                                                          'disruption'])

                aa_l, aa_ins, aa_dels, aa_disr = [ i if i else '' for i in [aa, list_ins, list_dels, disr] ]
                
                if genes:
                    g_hom[k][field] = aa_l
                    g_hom[k][field2] = len(aa_l)
                    g_hom[k][field3] = ';'.join(aa_ins)
                    g_hom[k][field4] = ';'.join(aa_dels)
                    g_hom[k][field5] = aa_disr
                    
                else:
                    d_hom[k][field] = aa_l 
                    d_hom[k][field2] = len(aa_l)
                    d_hom[k][field3] = ';'.join(aa_ins)
                    d_hom[k][field4] = ';'.join(aa_dels)
                    d_hom[k][field5] = aa_disr
                


    return fasta_IGHD_list, junction_seqs, d_hom, g_hom


def subsetFunction(fhand, ref, in_folder, out, gene=False):

    """"""
        
    bam_out = os.path.join(out, (fhand.replace('-sorted.bam', '_' + \
                                               ref.replace('-D','D-') +
                                               '-sorted.bam')))
    bam = os.path.join(in_folder, fhand)

    if gene:
        bam_out_gene = os.path.join(out, (fhand.replace('-sorted.bam', '_' + \
                                                        ref.replace('-D','D-') +
                                                        '-gene-sorted.bam')))
        
        subset = ("samtools view -H {0} > {2}; samtools view {0} -F4"
                  "| awk '{{ if ($3 ~ /{1}D*-*/) {{print}} }}' "
                  ">> {2}; samtools view {2} -b -o {3}; samtools index {3};"
                  " rm {2}").format(bam, ref, bam_out, bam_out_gene)
        
        bam_out = bam_out_gene
        
    else:
        subset = ('samtools view -h {0} -F4 {1} -b -o {2}; samtools'
                  ' index {2}').format(bam, ref.replace('-D','D-'), bam_out)

    return bam_out, subset


def mergeFunction(bam, ref, out, outV, outJ):

    """"""
    out_name = '{}_{}-merged-sorted.bam'.format(bam.replace('-sorted.bam',''),
                                                ref)

    out_merged_path = os.path.join(out, out_name)

    CMD = ('samtools merge {0} {1} {2} -f; '
           'samtools index {0}').format(out_merged_path, outV,
                                        outJ)
    
    return out_merged_path, CMD


def subset_merge(d, g, CMD_subset, CMD_merge, out, out_merged,
                 bamsV, bamsJ):

    """"""
    l_merged, j, v = ([], [], [])
    fCMD_subset = open_file(CMD_subset, mode='add')
    fCMD_merge = open_file(CMD_merge, mode='add')

    for e in d:
        
        vr, vgene, jr, vjr, vjr_gene = ([], [], [], [], [])
        sample = e.split('_')[0] #clonal
        refV = e.split('_')[-1].replace('D-','-D') #clonal
        geneV = '-'.join(refV.replace('D-','-D').split('-')[:-1])
        refVJ = d[e]['IGHV-J'].strip()
        refJ = d[e]['J_assigned'].strip()
        VJgene = '{}-gene_{}'.format(geneV, refJ)
        bamname = '{}-sorted.bam'.format(sample)
        
        # subset IGHV
        if refV != '' and refV not in vr:
            vr.append(refV)

            ## BAM subset
            bam_outV, subsetV = subsetFunction(bamname,
                                               refV, bamsV, out)
            
            if bam_outV not in v:
                v.append(bam_outV)
                fCMD_subset.write(subsetV + '\n')
                
        if e in g:
            
            # subset IGHV gene
            if geneV != '' and geneV not in vgene:
                vgene.append(geneV)

                ## BAM subset
                bam_outgV, subsetVgene = subsetFunction(bamname,
                                                        geneV, bamsV,
                                                        out, gene=True)
                if bam_outgV not in v:
                    v.append(bam_outgV)
                    fCMD_subset.write(subsetVgene + '\n')


        # subset IGHJ
        if refJ != '' and refJ not in jr:
            jr.append(refJ)

            ## BAM subset
            bam_outJ, subsetJ = subsetFunction(bamname,
                                               refJ, bamsJ,
                                               out)
            
            if bam_outJ not in j:
                j.append(bam_outJ)
                fCMD_subset.write(subsetJ + '\n')

        # BAM merge of IGHV and IGHJ reads
        if refV != '' and refJ != '' and refVJ not in vjr:
            vjr.append(refVJ)

            ## BAM merge
            out_merged_path, merge = mergeFunction(bamname, refVJ,
                                                   out_merged, bam_outV,
                                                   bam_outJ)
            l_merged.append(out_merged_path)
            fCMD_merge.write(merge + '\n')
            
        if e in g:
            if geneV != '' and refJ != '' and VJgene not in vjr_gene:
                vjr_gene.append(VJgene)

                ## BAM merge
                out_merged_path, merge = mergeFunction(bamname, VJgene,
                                                       out_merged, bam_outgV,
                                                       bam_outJ)
                
                l_merged.append(out_merged_path)
                fCMD_merge.write(merge + '\n')
                
    return l_merged


def merged_list(d, out_merged):

    """"""
    l_merged = []
    for e in d:
        sample, refV = e.split('_')
        refVJ = d[e]['IGHV-J'].strip()
        refJ = d[e]['J_assigned'].strip()
        bamname = sample + '-sorted.bam'
        if refV != '' and refJ != '':
            out_name = '{}_{}-merged-sorted.bam'.format(bamname.replace('-sorted.bam',''),
                                                        refVJ)
            out_merged_path = os.path.join(out_merged, out_name)
            l_merged.append(out_merged_path)

    return l_merged


def mapping_IGHVJ_assigned(d, g, bamsV, bamsJ, out_subset, out_merged):

    """"""
    overwrite = do_overwrite(out_subset, '.bam')
    if overwrite:
        log.debug('Starting BAM subset for IGHV and IGHJ represented genes.' 
                  ' Out folder %s', out_subset)
        path_CMD_subset, path_LOG_subset = create_CMD('bamsubset', out_subset)
        path_CMD_merge, path_LOG_merge = create_CMD('bamsmerge', out_subset)

        list_merged = subset_merge(d, g, path_CMD_subset, path_CMD_merge, 
                                   out_subset, out_merged, bamsV, bamsJ)
        parallel(path_CMD_subset, proc, path_LOG_subset)
        parallel(path_CMD_merge, proc, path_LOG_merge)

    else:
        log.warning('There are BAM files inside %s. Use mode overwrite if you '
                    'want to create BAM files in this folder', out_merged)
        list_merged = merged_list(d, out_merged)

    return list_merged


def consensus_sequence_annotation(dictionary, g, fasta_path, tag):

    """"""
    if clonal:
        read = os.path.basename(fasta_path).split('_')[1]
        sample = '{}_{}'.format(os.path.basename(fasta_path).split('_')[0], read) 
        Vgene = os.path.basename(fasta_path).split('_')[2].replace('-fb.fa','')
        
    else:
        sample = os.path.basename(fasta_path).split('_')[0]
        Vgene = os.path.basename(fasta_path).split('_')[1]
    
    if 'gene' in Vgene:
        for a in g:
            if sample == a.split('_')[0]:
                if (Vgene.replace('gene','')) in a.replace('D-','-D'):
                    k = a
    else:
        k = '{}_{}'.format(sample, Vgene)
        if clonal:
            
            k = '{}_{}'.format(sample.replace('123456789', '-'), Vgene)
            
    seq = read_file_simply(fasta_path)
    sequence = ''.join(seq[1:])
    dictionary[k][tag] = len(sequence)
    if k in g:
        g[k][tag] = len(sequence)

    return dictionary
 

def prepare_complete_consensus(tag, list_bams, out_folder, folder_vcfs, 
                               folder_rearrangement_references):
    
    """"""
    # iterate in d_hom
    overwrite = do_overwrite(out_folder, tag + '.fa')
    consensus_sequences = []
    vcfs = []
    if overwrite:
        log.debug('Starting complete rearrangement consensus sequence. '
                  'Out folder is %s', out_folder)
        path_CMD_consensus, path_LOG_consensus = create_CMD('complete_consensus', out_folder)
        fcmd = open_file(path_CMD_consensus, mode='write')

        for bam in list_bams:
            out_name = os.path.basename(bam).replace('-sorted.bam','')
            sample_name = out_name.split('_')[0]
            current_ref = '_'.join(out_name.split('_')[1:]).replace('-noFR3', '')
            complete_ref = '_'.join(out_name.split('_')[1:])
            consensus_path_fb, vcf_path = consensus_sequence_complete(current_ref, complete_ref, out_folder, 
                                                                      folder_rearrangement_references, bam, 
                                                                      sample_name, fcmd, folder_vcfs)

            # add consensus sequences paths to a list
            consensus_sequences.append(consensus_path_fb)
            vcfs.append(vcf_path)

        fcmd.close()
        log.info('Executing %s; generating complete consensus sequences for each '
                 'sample and probable allele', path_CMD_consensus)
        parallel(path_CMD_consensus, proc, path_LOG_consensus)
        parallel_resume(path_CMD_consensus, proc, path_LOG_consensus)

    else:
        log.warning('There are already consensus sequences files in %s. '
                    'Use mode overwrite if you want to overwrite the existing file.', out_folder)
        consensus_sequences = glob.glob(str(out_folder) + ('/*' + tag + '*.fa'))
        vcfs = glob.glob(str(folder_vcfs) + (tag + '/*' + tag + '*.vcf.gz'))


    return consensus_sequences, vcfs


def prepare_consensus_sequence_annotation(dictionary, g, list_consensus):

    """annotate consensus sequence length to d_hom"""
    for path in list_consensus:

        if path.endswith('-fb.fa'):
            tag = 'consensus_length'

            consensus_sequence_annotation(dictionary, g, path, tag)


def bam2fastq(fcmd, bam_name, out):

    """"""
    fastq_name = os.path.basename(bam_name).replace('-merged-sorted.bam', '.fastq.gz')
    fastq = os.path.join(out, fastq_name)
    CMD = 'bamtools convert -format fastq -in {} | gzip > {}'.format(bam_name, fastq)
    fcmd.write(CMD + '\n')

    return fastq


def merged_fastqs_list(l, out):

    """"""
    out_l = []
    for name in l:
        fastq = os.path.basename(name).replace('-merged-sorted.bam', '.fastq.gz')
        fastq_path = os.path.join(out, fastq)
        out_l.append(fastq_path)

    return out_l


def merged2fastq(l, out):

    """"""
    overwrite = do_overwrite(out, '.fastq.gz')
    if overwrite:
        list_fastqs = []
        log.debug('Starting BAM subset for IGHV and IGHJ merged read conversion to FASTQ. '
                  'Out folder %s', out)

        path_CMD_fastq, path_LOG_fastq = create_CMD('merged2fastq', out)
        fcmd = open_file(path_CMD_fastq, mode='write')
        for bam in l:
            fastq_name = bam2fastq(fcmd, bam, out)
            list_fastqs.append(fastq_name)
        fcmd.close()
        parallel(path_CMD_fastq, proc, path_LOG_fastq)
    else:
        log.warning('There are FASTQ files inside %s. Use mode overwrite if you want '
                    'to create FASTQ files in this folder', out)
        list_fastqs = merged_fastqs_list(l, out)

    return list_fastqs


def rearrangement_bams_list(l, out):

    """"""
    new_l = []
    for fastq in l:
        bam_name = os.path.basename(fastq).replace('.fastq.gz','-sorted.bam')
        bam_path = os.path.join(out, bam_name)
        new_l.append(bam_path)

    return new_l


def specific_rearrangement_bwamem(fastq, out, fcmd, ref_dict, out_ref, g):

    """"""
    sample_name = os.path.basename(fastq).split('_')[0]
    ref_name = '_'.join(os.path.basename(fastq).replace('.fastq.gz','').split('_')[1:])
    if 'gene' in ref_name:
        for s in g:
            if sample_name == s.split('_')[0]:
                if ref_name.split('_')[0].replace('gene', '') in s.replace('D-','-D'):
                    allele = s.split('_')[1] + '_' + ref_name.split('_')[1]
                    
    else:
        allele = ref_name
        
    seq = ref_dict[allele]
    bam_name = os.path.basename(fastq).replace('.fastq.gz','-sorted.bam')
    bam_path = os.path.join(out, bam_name)
    ref_path = os.path.join(out_ref, (ref_name + '.fa'))
    make_reference(ref_path, ref_name, seq, mode='write')
    CMD = ('bwa index {1}; bwa mem -t4-M -L 50 -R '
           '\"@RG\\tID:{0}\\tPL:ILLUMINA\\tSM:{0}\\tDS:ref={1}\\tCN:UGDG\\tDT:{4}\\tPU:{0}\" '
           '{1} -p {2} | samtools view -b - | samtools sort -o {3} -; samtools index {3}').format(sample_name, 
                                                                                                  ref_path, fastq, 
                                                                                                  bam_path, cdate)

    fcmd.write(CMD + '\n')

    return bam_path


def uniq_fastq_script(fastq, fcmd):

    """"""
    script = os.path.join(path_IGH_scripts, 'dedup.py')
    CMD = 'python3.5 {} {}'.format(script, fastq)
    fcmd.write(CMD + '\n')
    

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
            if m[read].split('/')[1] == '1':
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


def remove_FR3(BAM, f):

    """"""
    out_BAM = BAM.replace('-sorted.bam', '-noFR3-sorted.bam')
    ref = '_'.join(os.path.basename(BAM).split('_')[1:]).replace('-sorted.bam', '')
    CMD = 'samtools view {0} -F16 {1}:1-200 -b -o {2}; samtools index {2}'.format(BAM, ref, out_BAM)
    f.write(CMD + '\n')

    return out_BAM
    
    
def specific_rearrangement_mapping(l, out, ref_dict, out_ref, g):

    """"""
    overwrite = do_overwrite(out, '.bam')
    if overwrite:
        list_bams = []
        list_bams_noFR3 = []
        log.debug('Starting mapping for specific rearrangements. Out folder %s', out)
        path_CMD_mapping, path_LOG_mapping = create_CMD('specific-mapping', out)
        path_CMD_removeFR3, path_LOG_removeFR3 = create_CMD('remove-FR3', out)
        path_CMD_uniq, path_LOG_uniq = create_CMD('uniq-fastq', out)
        fcmd0 = open_file(path_CMD_uniq, mode='write')
        fcmd = open_file(path_CMD_mapping, mode='write')
        fcmd2 = open_file(path_CMD_removeFR3, mode='write')
        for fastq in l:
            
            uniq_fastq_script(fastq, fcmd0)
            bam_name = specific_rearrangement_bwamem(fastq, out, fcmd, ref_dict, out_ref, g)
            list_bams.append(bam_name)

            ## create BAM without FR3
            bam_noFR3_name = remove_FR3(bam_name, fcmd2)
            list_bams_noFR3.append(bam_noFR3_name)

        fcmd0.close()
        fcmd.close()
        fcmd2.close()
        parallel(path_CMD_uniq, proc, path_LOG_uniq)
        
        if len(list_bams) >= 4:
            j = int(int(proc)/4)
        else:
            j = proc
            
        parallel(path_CMD_mapping, j, path_LOG_mapping)
        parallel(path_CMD_removeFR3, proc, path_LOG_removeFR3)
        
    else:
        log.warning('There are BAM files inside %s. Use mode overwrite if you '
                    'want to create BAM files in this folder', out)
        lbams = glob.glob(str(out) + '/*-sorted.bam')
        list_bams, list_bams_noFR3 = [], []
        for l in lbams:
            if 'noFR3' in l:
                list_bams_noFR3.append(l)
            else:
                list_bams.append(l)

    return list_bams, list_bams_noFR3


def sequences_overlap(s1, s2):

    """"""
    d = difflib.SequenceMatcher(None, s1, s2)
    match = max(d.get_matching_blocks(), key=lambda x:x[2])
    i, j, k = match
    #i: pos where match starts in s1
    #j: pos where match starts in s1
    #k: number of characters overlapping

    return i, j, k


def nspscalculation(junction, vseq, jseq, dseq):

    """"""
    
    junction_start, d_start, length_juncd = sequences_overlap(junction, dseq)
    vreg_seq = junction[:sequences_overlap(junction, dseq)[0]]
    jreg_seq = junction[(sequences_overlap(junction, dseq)[0] + sequences_overlap(junction, dseq)[2]):]

    n = junction[:sequences_overlap(vreg_seq, vseq)[2]]
    p = junction[:-(sequences_overlap(jreg_seq, jseq)[2])][-(sequences_overlap(jreg_seq, jseq)[0]):]

    return (n, p)


def check_VH(seq, bam_name, join_info):

    """"""
    sample = os.path.basename(bam_name).split('_')[0]
    Vallele = os.path.basename(bam_name).split('_')[1]
    ftype = os.path.join(join_info, ('info_bams_{}.txt'.format(sample)))
    CMD = 'grep -A7 "{}" {} | grep {} | wc -l'.format(seq, ftype, Vallele)
    CMD2 = 'grep -A7 "{}" {} | grep "IGHV" | wc -l'.format(seq, ftype)
    proc = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE)
    n = float(proc.communicate()[0].strip())
    proc2 = subprocess.Popen(CMD2, shell=True, stdout=subprocess.PIPE)
    n2 = float(proc2.communicate()[0].strip())
    
    if n2 > 0:
        if n/n2*100 > 1:

            return True
        else:
            return False
    

def finding_new_junction(file_reads, vcf_complete_path, bam, refV_seq, info_folder):

    """Store all unique reads in a file and start looking for CDR3"""
    CMD3 = ('samtools view {0} -F4 | cut -f10 | sort | uniq -c | sort -n > {1};'
            'cat {1} | tail -n 1').format(bam, file_reads)
    totalCMD = 'samtools view {0} -F4 | wc -l'.format(bam)
    proc1 = subprocess.Popen(totalCMD, shell=True, stdout=subprocess.PIPE)
    totalreads = float(proc1.communicate()[0].strip())
    proce = subprocess.Popen(CMD3, shell=True, stdout=subprocess.PIPE)
    ns = proce.communicate()[0].strip().split()
    ## if BAM file is empty we do not do this
    if len(ns) > 1:
        major_readn, major_read = ns
        # we define major_read as the sequence with more reads first
        list_majors = read_file_simply(file_reads)
        l_junction = None
        n = iter((e.strip().split() for e in list_majors[::-1]))
        try:
            while not l_junction:
                ## will return an empty junction if CDR3 not found
                l_next = next(n)
                
                ## include sequence and sample name to search VH info
                if float(l_next[0]) >= (totalreads/100)*3:
                    if check_VH(l_next[1], bam, info_folder):
                        IGHD, pos_IGHD, list_insertions, list_deletions, disruption = parse_vcf_IGHD(vcf_complete_path,
                                                                                                 l_next[1], 290)
                        l_junction, start, end, prod = consensus2CDR3.cdr3_extraction(l_next[1], mincys=3)
                
                        major_readn, major_read = l_next
                        if IGHD == '' or not IGHD:
                            IGHD = calculateD(major_read, start, end)
                    else:
                        (l_junction, IGHD, start, end,
                         major_read, major_readn, prod) = [None] * 7
                else:
                    (l_junction, IGHD, start, end,
                     major_read, major_readn, prod) = [None] * 7
                
        except StopIteration:
            pass
        
        finally:
            del n
    else:
        l_junction, IGHD, start, end, major_read, major_readn, prod = [None] * 7

    
    return l_junction, IGHD, start, end, major_read, major_readn, prod
    

def major_reads_junction(bam, seq, file_reads):

    """Representation of the original junction"""
    percent_junction = 0
    ## CMDs
    CMD = 'samtools view {} -F4 | grep -i {} | wc -l'.format(bam, seq)
    CMD2 = 'samtools view {} -F4 | wc -l'.format(bam)

    reads_junction = process_communicate(CMD)[0]
    
    # if the junction seq is represented we calculate in which percentage
    if int(reads_junction) > 0:
        reads_mapped = process_communicate(CMD2)[0]
        percent_junction = float(reads_junction)/float(reads_mapped)
                                                                 
    return percent_junction


def calculateD(clean_seq, real_start, real_end):

    """"""
    if all([clean_seq, real_start, real_end]):
        seq = clean_seq[real_start : real_end]
    else:
        seq = ''
        
    return seq


def bams_specific_read(list_bams, out, vcf_folder,
                       homology_folder, Vdict, Jdict, d, info_folder):

    """"""
    overwrite = do_overwrite(out, '.txt')
    if overwrite:
        list_seqs = []
        log.debug('Starting looking for unique sequences for specific rearrangements.'
                  ' Out folder %s', out)
        
        fasta_IGHD_list = []
        fasta_list = []
        def_list_bams = [b for b in list_bams if not 'gene' in b]
        for BAM in def_list_bams:
            
            fastaname = os.path.basename(BAM).replace('.bam', '-IGHD.fasta')
            fastaD_path = os.path.join(out, fastaname)
            fasta_path = fastaD_path.replace('-IGHD', '')
            knam, refv, rest = os.path.basename(BAM).split('_')
            kname = '{}_{}'.format(knam, refv)
            vcf_incomplete_path = os.path.join(vcf_folder, kname)
            vcf_complete_path = glob.glob(vcf_incomplete_path + '*.vcf.gz')[0]
            refV_seq = Vdict[refv]
            refJ = rest.replace('-sorted.bam','')
            
            refJ_seq = Jdict[refJ]
            
            #seq_junction = junctiond[kname]
            
            f_reads = os.path.join(out, (kname + '_uniqueread_counts.txt'))
            list_seqs.append(f_reads)
            new_junction, IGHD, start, end, seq, seqn, prod = finding_new_junction(f_reads,
                                                                             vcf_complete_path,
                                                                                   BAM, refV_seq,
                                                                                   info_folder)
            if not seq:
                read_CDR3 = ''
            else:
                read_CDR3 = seq[start:end]
            
            if new_junction:
                major_read_np = nspscalculation(new_junction, refV_seq,
                                                refJ_seq, IGHD)
            ## define percent junction for non junction rearrangement
            ## LOOK FOR THE JUNCTION IN BAM ONLY IF IT EXISTS
            #if isinstance(seq_junction, list):
            #    list_junction = [t for t in seq_junction if t != '']
                
            #    for s in list_junction:
                    
            #        original_junction = s
            #        # if the original region is > 20% in BAMS
            #        percent_junction = major_reads_junction(BAM, s, f_reads)
            #        # ori_read_np = nspscalculation(seq, refV_seq, refJ_seq, IGHD)
            #        # put in IGHD extraction part
                    
            #elif seq_junction != '':
            #    original_junction = seq_junction
            #    # if the original region is > 20% in BAMS
            #    percent_junction = major_reads_junction(BAM, seq_junction, f_reads)
            #    # ori_read_np = nspscalculation(seq, refV_seq, refJ_seq, IGHD)
                
                

            #else:
            #    original_junction = ''
            #    percent_junction = 0
               
            #d[kname]['ori_junction'] = original_junction
            #d[kname]['perc_ori_junction'] = percent_junction
            if not IGHD:
                IGHD = ''
            d[kname]['new_seq'] = seq
            d[kname]['nreads_new_seq'] = seqn
            d[kname]['new_CDR3'] = new_junction
            d[kname]['new_IGHD'] = IGHD
            CDR32fasta([IGHD], out, fastaD_path)
            CDR32fasta([read_CDR3], out, fasta_path)
            fasta_IGHD_list.append(fastaD_path)
            fasta_list.append(fasta_path)
            
            d[kname]['prod'] = prod
    else:
        list_seqs =  glob.glob(str(out) + '/*_uniqueread_counts.txt')
        fasta_IGHD_list = glob.glob(str(out) + '/*-D.fa')


    return list_seqs, d, fasta_IGHD_list

    # write:     
    # - original_junction: junction sequence is there is one
    # - percent_junction: percent of reads original_junction
    # - seq: most abundant sequence with a productive junction
    # - seqn: number of reads seq
    # - new_junction: junction found in seq
    # - IGHD: IGHD in the new junction
    # - prod: if new_junction is productive 
    
    
def CDR3_short_blast(out_blast, fastq, sample_name, fcmd):

    """"""
    # paths and commands
    ## blast files
    ## check first if IGHD is empty
    CMD2 = 'wc -l {}'.format(fastq)
    proce = subprocess.Popen(CMD2, shell=True, stdout=subprocess.PIPE)
    out = proce.communicate()[0].decode("utf-8").split(' ')[0]
    # if IGHD fasta is empty, we use complete CDR3 sequence
    if out == '0':
        fastq = fastq.replace('-IGHD', '')
    # with sorting we are ordering BLAST file by bitscore
    CMD = ('blastn -db IMGT_D.fa -task blastn-short -evalue 1000 -dust no -strand plus -word_size 5 ' 
           '-query {} -outfmt 6 |  sort -k12 -nr > {}').format(fastq, out_blast)
    fcmd.write(CMD + '\n')


def CDR3_emboss(out_emboss, fastq, fcmd, refD, emboss_files):

    """Emboss water alignment for each IGHD fasta with each IGHD allele"""
    allelesD = VJref_dictionary(refD)
    for k in allelesD:
        nameD = k.split('|')[1]
        out_file = '{}_{}.water'.format(out_emboss, nameD)
        ## check first if IGHD is empty
        CMD2 = 'wc -l {}'.format(fastq)
        proce = subprocess.Popen(CMD2, shell=True, stdout=subprocess.PIPE)
        out = proce.communicate()[0].decode("utf-8").split(' ')[0]
        # if IGHD fasta is empty, we use complete CDR3 sequence
        if out == '0':
            fastq = fastq.replace('-IGHD', '')
            CMD2 = 'wc -l {}'.format(fastq)
            proce = subprocess.Popen(CMD2, shell=True, stdout=subprocess.PIPE)
            out = proce.communicate()[0].decode("utf-8").split(' ')[0]
            # if both IGHD and CDR3 fasta are empty, the loop for that rearrangement will be interrupted
            if out == '0':
                break
            else:
                emboss_files.append(out_file)
        else:
            emboss_files.append(out_file)

        CMD = ('emboss_water -asequence {} -bsequence <(echo -e ">{}\\n{}") -gapopen 10.0 '
               '-gapextend 0.5 -outfile {} -aformat pair').format(fastq, nameD, allelesD[k], out_file)

        fcmd.write(CMD + '\n')

    return emboss_files


def parse_alignment_CDR3_emboss(fasta_CDR3_list, out_folder, emboss_files):

    """"""
    for fastq in fasta_CDR3_list:
        # open resume file and parse alignments
        sample_name = os.path.basename(fastq).replace('.fasta', '')
        out_resume_emboss = os.path.basename(fastq).replace('.fasta', '-EW')
        out_emboss = os.path.join(out_folder, out_resume_emboss)
        unsorted = out_emboss.replace('-EW','-EWunsorted')
        fresume = open_file(unsorted, mode='write')

        for e in emboss_files:
            eb = os.path.basename(e)
            ebpre = '_'.join(eb.split('_')[:-1]).replace('-EW','')
            if ebpre == sample_name:
                IGHD = eb.split('_')[-1].replace('.water', '')
                identity, length, score = parse_emboss(e)
                fresume.write(IGHD + ',' + identity + ',' + 
                              length.split(' ')[2].replace(";","") + 
                              ',' + score + '\n')

        fresume.close()
        # Order the sequences per identity and length
        CMD = 'sort -k4 -t"," -nr {} > {}'.format(unsorted, out_emboss)
        execute(CMD)


def prepare_CDR3_emboss(fasta_CDR3_list, out_folder, refD):

    """"""
    emboss_files = []
    EW_list = []
    overwrite = do_overwrite(out_folder, 'IGHD-EW')
    if overwrite:
        log.debug('Starting CDR3 emboss water alignment. Out folder %s', out_folder)
        path_CMD_emboss, path_LOG_emboss = create_CMD('CDR3_emboss', out_folder)
        fcmd = open_file(path_CMD_emboss, mode='write')

        ## generate CMD and align each sample against all D alleles
        for fastq in fasta_CDR3_list:
            # prepare emboss files
            #sample_name = os.path.basename(fastq).split('_')[0]
            #IGHV = os.path.basename(fastq).split('_')[1]
            out_name_emboss = os.path.basename(fastq).replace('.fasta', '-EW')
            out_emboss = os.path.join(out_folder, out_name_emboss)
            EW_list.append(out_emboss)
            emboss_files = CDR3_emboss(out_emboss, fastq, fcmd, refD, emboss_files)

        fcmd.close()
        parallel(path_CMD_emboss, proc, path_LOG_emboss)

        ## parse alignments per sample
        parse_alignment_CDR3_emboss(fasta_CDR3_list, out_folder, emboss_files)

    else:
        log.warning('There are BLAST files inside %s. Use mode overwrite if you want '
                    'to create BLAST files in this folder', out_folder)
        emboss_files = glob.glob(str(out_folder) + '/*EW_IGHD*')
        EW_list = glob.glob(str(out_folder) + '/*IGHD-EW')

    
    return emboss_files, EW_list


def prepare_CDR3_blast(fasta_CDR3_list, out_folder):

    """short BLAST to identify IGHD"""
    blast_files = []
    overwrite = do_overwrite(out_folder, 'IGHD')
    if overwrite:
        log.debug('Starting CDR3 blast. Out folder %s', out_folder)
        path_CMD_blast, path_LOG_blast = create_CMD('CDR3_blast', out_folder)
        fcmd = open_file(path_CMD_blast, mode='write')

        for fastq in fasta_CDR3_list:
            # prepare blast files
            sample_name = os.path.basename(fastq).split('_')[0]
            IGHV = os.path.basename(fastq).split('_')[1]
            out_name = os.path.basename(fastq).replace('.fasta', '')
            out_blast = os.path.join(out_folder, out_name)
            CDR3_short_blast(out_blast, fastq, sample_name, fcmd)

            ## add blast tabular results file to list
            blast_files.append(out_blast)

        fcmd.close()
        log.info('Executing %s; generating results tabular format files from '
                 'EMBOSS water performance with CDR3 sequences extracted for '
                 'each sample and probable allele', path_CMD_blast)

        parallel(path_CMD_blast, proc, path_LOG_blast)

    else:
        log.warning('There are water files inside %s. Use mode overwrite if you '
                    'want to create BLAST files in this folder', out_folder)
        blast_files = glob.glob(str(out_folder) + '/*IGHD')
        

    return blast_files


def emboss_annotation(emboss_list, d_hom, g_hom):

    """"""
    field = ''
    for emboss in emboss_list:
        if clonal:
            sample_name = '_'.join(os.path.basename(emboss).split('_')[0:2])
            sample_name = sample_name.replace('123456789','-')
            IGHV = os.path.basename(emboss).split('_')[2].replace('-IGHD-EW','')
        else:
            sample_name = os.path.basename(emboss).split('_')[0]
            IGHV = os.path.basename(emboss).split('_')[1]
            print('we')
            print('sample_name', 'IGHV')
        if 'gene' in IGHV:
            for a in g_hom:
                if sample_name == a.split('_')[0]:
                    if (IGHV.replace('gene','')) in a:
                        k = a
        else:
            k = sample_name + '_' + IGHV

        field = 'IGHD_emboss'

        fblast = read_file_simply(emboss)
        
        if fblast != []:
            if len(fblast) >= 3:
                CDR3 = ' '.join(fblast[0:3]).replace(',',' ')
            elif len(fblast) >=2:
                CDR3 = ' '.join(fblast[0:2]).replace(',',' ')
            elif len(fblast) >= 1:
                CDR3 = fblast[0].replace(',',' ')
        else:
            CDR3 = ''

        if len(IGHV.split('-')) == 2:
            if k in g_hom:
                g_hom[k][field] = CDR3
        else:
            d_hom[k][field] = CDR3


    return d_hom, g_hom


def blast_annotation(blast_list, d_hom, g_hom):

    """"""
    for blast in blast_list:
        if clonal:
            sample_name = '_'.join(os.path.basename(blast).split('_')[0:2])
            sample_name = sample_name.replace('123456789','-')
            IGHV = os.path.basename(blast).split('_')[2].replace('-fb-IGHD','')
        else:
            sample_name = os.path.basename(blast).split('_')[0]
            IGHV = os.path.basename(blast).split('_')[1]

        if 'gene' in IGHV:
            for a in g_hom:
                if sample_name == a.split('_')[0]:
                    if (IGHV.replace('gene','')) in a:
                        k = a
        else:
            k = sample_name + '_' + IGHV

        if 'fb' in os.path.basename(blast):
            field = 'IGHD'

        fblast = read_file_simply(blast)
        if fblast != []:
            if len(fblast) >= 3:
                CDR3 = (fblast[0].split('\t')[1].split('|')[1] + ' ' + 
                        ' '.join(fblast[0].split('\t')[2:4]) +
                        ';' + fblast[1].split('\t')[1].split('|')[1] + 
                        ' ' + ' '.join(fblast[1].split('\t')[2:4]) +
                        ';' + fblast[2].split('\t')[1].split('|')[1] + 
                        ' ' + ' '.join(fblast[2].split('\t')[2:4]))
            elif len(fblast) >=2:
                CDR3 = (fblast[0].split('\t')[1].split('|')[1] + ' ' + 
                        ' '.join(fblast[0].split('\t')[2:4]) +
                        ';' + fblast[1].split('\t')[1].split('|')[1] + 
                        ' ' + ' '.join(fblast[1].split('\t')[2:4]))
            elif len(fblast) >= 1:
                CDR3 = (fblast[0].split('\t')[1].split('|')[1] + ' ' + 
                        ' '.join(fblast[0].split('\t')[2:4]))
        else:
            CDR3 = ''

        if len(IGHV.split('-')) == 2:
            if k in g_hom:
                g_hom[k][field] = CDR3
        else:
            d_hom[k][field] = CDR3


    return d_hom, g_hom


def main():

    """Main script for pipeline intro"""
    # Get arguments
    args = arguments_parser_IGH.parse_args()

    ## Define logging
    log_setup(args)

    # define global variables
    global cdate
    cdate = time.strftime("%d-%m-%Y_%H-%M-%S")
    global proc
    global owrt
    global clonal
    global frec
    global minalt
    global basal
    global mincoverage
    global TMPDIR
    out_folder = os.path.abspath(args.dest)
    TMPDIR = out_folder
    # analysis mode
    if args.pipeline:
        fastqs_folder = os.path.abspath(args.inputf)
        owrt = args.owrt
        proc = args.processes
        # Define folder exclusion from FASTQ seeking
        if args.exclude:
            exclud = args.exclude
        else:
            exclud = 'None'
        # basal filter
        if args.basal:
            basal = 3.6
        else:
            basal = 0
        # clonality analysis?
        if args.clonality:
            clonal = True
            basal = 0
        else:
            clonal = False

        # kind of analysis (Sanger or NGS)
        if args.sanger:
            frec = '0.2'
            minalt = '1'
            mincoverage = '0'
        else:
            frec = '0.5'
            minalt = '2'
            mincoverage = '50'

    # test option
    if args.test:
        fastqs_folder = path_test #"/nfs/hpcugdga/IGHS-motivos_20171205/fastq_test"
        owrt = args.owrt
        proc = '10'
        exclud = 'None'
        basal = 3.6
        frec = '0.5'
        minalt = '2'
        mincoverage = '50'
        # clonality analysis?
        if args.clonality:
            clonal = True
            basal = 0
        else:
            clonal = False


    # Merging whole process
    (fastq_dictionary, merged_folder, pairs_raw,
     pairs_trimming, trimming_io, extensions) = fastq_merge.whole_merging_process(fastqs_folder, 
                                                                                  out_folder, proc, cdate, exclud)
    ext = set(extensions)

    # Trimming
    ## Perform quality trimming plus primers or only quality trimming
    trimmed_folder = os.path.join(out_folder, 'trimmed')
    create_dir(trimmed_folder)

    if args.primers:
        adapters = args.primers
        trimmed5_folder = os.path.join(trimmed_folder, 'trimmed5')
        
        if args.sides:
            create_dir(trimmed5_folder)
            prepare_trimming(trimmed5_folder, fastq_dictionary, merged_folder,fastqs_folder, 
                             ext, pairs_raw, adapters, args.sides)
            prepare_trimming(trimmed_folder, fastq_dictionary, trimmed5_folder, fastqs_folder, 
                             ext, pairs_raw, adapters, args.sides, prog='bbduk', mode='right')
        else:
            prepare_trimming(trimmed_folder, fastq_dictionary, merged_folder,fastqs_folder, 
                             ext, pairs_raw, adapters, args.sides)
    else:
        adapters = 0
        prepare_trimming(trimmed_folder, fastq_dictionary, merged_folder, fastqs_folder, 
                         ext, pairs_raw, adapters, args.sides, prog='seqtk')


    # Merge R1 and R2 using FLASH (only for clonality mode)
    if clonal:

        fastq_merged_folder = os.path.join(out_folder, 'trimmed_merged_clonality')
        create_dir(fastq_merged_folder)
        new_pairs_flash, flash_boolean = prepare_flash(fastq_merged_folder, 
                                                       trimmed_folder, pairs_trimming)
        if flash_boolean:
            trimmed_folder = fastq_merged_folder
        else:
            trimmed_folder = trimmed_folder
        pairs_trimming = new_pairs_flash

    # Alignment with BWA mem (autodetect paired or single reads)
    ## V bams
    alignment_folderV = os.path.join(out_folder, 'bamsV')
    create_dir(alignment_folderV)
    refV = os.path.join(IGH_references, 'IMGT_V-oneline-NOstars-simple.fa')
    bam_listV = prepare_alignment(alignment_folderV, trimmed_folder, pairs_trimming, refV)
    ## D bams
    alignment_folderD = os.path.join(out_folder, 'bamsD')
    create_dir(alignment_folderD)
    refD = os.path.join(IGH_references, 'IMGT_D-nostars.fa')
    #bam_listD = prepare_alignment(alignment_folderD, trimmed_folder, pairs_trimming, refD)
    ## J bams
    alignment_folderJ = os.path.join(out_folder, 'bamsJ')
    create_dir(alignment_folderJ)
    refJ = os.path.join(IGH_references, 'IMGT_J-nostars-simple.fa')
    bam_listJ = prepare_alignment(alignment_folderJ, trimmed_folder, pairs_trimming, refJ)

    # Alignment statistics with samtools idxstats and probable regions calculation
    ## V bams
    results_folder = os.path.join(out_folder, 'results')
    create_dir(results_folder)
    stats_folderV = os.path.join(out_folder, 'statsV')
    create_dir(stats_folderV)
    statsV_list = prepare_stats(stats_folderV, fastq_dictionary, alignment_folderV, bam_listV)
    statsV_fof = list2file(statsV_list, 'statsV.fof', stats_folderV)
    path_probable_regionsV = probable_regions(statsV_fof, results_folder, 'probable_vregions.csv')

    ## J bams
    stats_folderJ = os.path.join(out_folder, 'statsJ')
    create_dir(stats_folderJ)
    statsJ_list = prepare_stats(stats_folderJ, fastq_dictionary, alignment_folderJ, bam_listJ)
    statsJ_fof = list2file(statsJ_list, 'statsJ.fof', stats_folderJ)
    path_probable_regionsJ = probable_regions(statsJ_fof, results_folder, 'probable_jregions.csv')

    # IGH BAM parsing
    ## We need a fof file of Vbams
    bamparse_folder = os.path.join(results_folder, 'bamparse')
    create_dir(bamparse_folder)
    info_folder = os.path.join(bamparse_folder, 'bamparsing_out')
    create_dir(info_folder)
    subtypes_folder = os.path.join(bamparse_folder, 'subtypes_resume')
    create_dir(subtypes_folder)
    consensus_complete = os.path.join(results_folder, 'consensus_complete')
    create_dir(consensus_complete)
    
    bamsV_fof = list2file(bam_listV, 'bamsV.fof', alignment_folderV)
    list_pars = bam_parsing(bamsV_fof, info_folder)
    ## Generate the subtype resume
    subtype_list = prepare_subtype(subtypes_folder, list_pars)
    # retrieve clonality region V and J files in clonal
    if clonal:
        clonality_V_file = os.path.join(info_folder, 'regions_clonalityV.csv')
        clonality_J_file = os.path.join(info_folder, 'regions_clonalityJ.csv')

    # paths
    ## Open output files
    f_homology = os.path.join(results_folder, 'homology_table.csv')
    f_hom = open_file(f_homology, mode='write')
    if not clonal:
        g_homology = os.path.join(results_folder, 'homology_table-Vgene.csv')
        m_hom = open_file(g_homology, mode='write')
    if args.basal:
        f_homology2 = os.path.join(results_folder, 'homology_table-filters.csv')
        f_hom2 = open_file(f_homology2, mode='write')

    # paths
    ## Open probable V regions and J regions resume
    regions_listV = read_file_simply(path_probable_regionsV)
    path_resume_regionsJ = path_probable_regionsJ.replace('probable_','resume_')
    regions_listJ = read_file_simply(path_resume_regionsJ)

    if clonal:
        regions_listV = read_file_simply(clonality_V_file)
        regions_listJ = read_file_simply(clonality_J_file)

    # header
    ## Write header to output
    if args.test:
        f_hom.write('Your result should be: MM,IGHV4-34-02_IGHJ6-02,'
                    'CARVITRASPGTDGRYGMDVW,IGHD3-16-02\n')

    ## write tab files header
    write_header(f_hom)
    
    if not clonal:
        write_header(m_hom)

    # create dictionary
    ## Create dictionary with V allele references
    ref_dict = VJref_dictionary(refV)
    ref_dictJ = VJref_dictionary(refJ)
    
    # paths
    ## Create folders for homology, consensus sequences for each V
    ## represented allele and CDR3 fastqs
    homology_folder = os.path.join(results_folder, 'homology')
    #homology_complete = os.path.join(results_folder, 'homology_complete')
    CDR3_fastas_folder = os.path.join(results_folder, 'CDR3_fastas')
    consensusV_folder = os.path.join(results_folder, 'consensus_Vregion')
    vcf_folder = os.path.join(out_folder, 'vcfs')
    alternativeJ_folder = os.path.join(out_folder, 'alternativeJ_assign')

    create_dir(homology_folder)
    create_dir(CDR3_fastas_folder)
    create_dir(consensusV_folder)
    create_dir(vcf_folder)
    create_dir(homology_folder)
    create_dir(alternativeJ_folder)

    # create a dictionary to store the values of the fields for each sample
    d_hom = {}
    g_hom = {}

    if not clonal and not args.sanger:
        # vregion filter
        (new_regions_listV, regions_listVgene, c_counts,
         c_counts_nofilter, nbasal) = filter_Vregion_basal(regions_listV,
                                                               bam_listV)
        new_regions_listJ = filter_Jregion_basal(regions_listJ, bam_listV)
            
    else:
        nbasal = basal
        

    # paths
    ## plots for clones proportion in samples
    plots_folder = os.path.join(results_folder, 'plots')
    plots_folder_nofilter = os.path.join(results_folder, 'nofilter_plots')
    create_dir(plots_folder)
    create_dir(plots_folder_nofilter)

    ### basal file
    basal_file = os.path.join(out_folder, 'basal_filter.txt')
    bas = open(basal_file, 'w')
    bas.write('The basal percent filtered is ' + str(nbasal) + '\n')
    bas.close()

    # pie plots V genes
    if args.basal and not clonal:
        pie_plots(c_counts, plots_folder)
    if not clonal and not args.sanger:
        pie_plots(c_counts_nofilter, plots_folder_nofilter)

    # d_hom creation
    d_hom = annotate_hom(bam_listV, regions_listJ, regions_listV, d_hom)
    if not clonal and not args.sanger:
        g_hom = annotate_hom(bam_listV, regions_listJ, regions_listVgene, g_hom)

    # Joined subtype annotation
    d_hom, g_hom = joined_annotation(d_hom, g_hom, subtype_list)
        # alternative J annotation (most cases of unique sanger sequences
        # do not work with the method above)
    if args.sanger:
        d_hom = alternative_J_assign(d_hom, pairs_trimming, alternativeJ_folder,
                                     trimmed_folder, refJ)

    # Mapping process against combinations of V and J
    ## new reference (create if it does not exist)
    references_folder = os.path.join(out_folder, 'references')
    create_dir(references_folder)
    refVJ = prepare_joined_references(references_folder, refV, refJ)

    # xpaths
    ## Mapping against IGHVJ-assigned
    subset_bams = os.path.join(out_folder, 'subset_bams')
    bams_merged = os.path.join(subset_bams, 'merged')
    folder_completevcfs = os.path.join(out_folder, 'vcfs_complete')
    merged2fastqs_folder = os.path.join(out_folder, 'merged2fastqs')
    specific_mapping = os.path.join(out_folder, 'bams_specific_rearrangements')
    rearrangement_refs = os.path.join(out_folder, 'rearrangement_refs')
    fblast_CDR3 = os.path.join(out_folder, 'CDR3_blast')
    uniq_bam_sequences = os.path.join(specific_mapping, 'uniq_sequences')
    
    # paths
    create_dir(subset_bams)
    create_dir(bams_merged)
    create_dir(folder_completevcfs)
    create_dir(merged2fastqs_folder)
    create_dir(specific_mapping)
    create_dir(rearrangement_refs)
    create_dir(fblast_CDR3)
    create_dir(uniq_bam_sequences)

    if clonal:
        # consensus must be a list with all fasta files per each sequence
        consensus_complete_list = glob.glob(str(consensus_complete) + '/*-fb.fa')

    else:
        list_merged = mapping_IGHVJ_assigned(d_hom, g_hom, alignment_folderV,
                                             alignment_folderJ, subset_bams, bams_merged)
        # Generate consensus read for each rearrangement
        ## Dict from VJ fasta reference file
        refVJ_dict = VJref_dictionary(refVJ)
        list_merged2fastq = merged2fastq(list_merged, merged2fastqs_folder)
        list_spec_bams, list_spec_noFR3_bams = specific_rearrangement_mapping(list_merged2fastq, specific_mapping,
                                                        refVJ_dict, rearrangement_refs, g_hom)
        
        # include step for vcf parsing for indels
        consensus_complete_list, vcfs_complete_list = prepare_complete_consensus('', list_spec_bams,
                                                                                 consensus_complete,
                                                                                 folder_completevcfs,
                                                                                 rearrangement_refs)

        vcf_parsing(results_folder, folder_completevcfs, vcfs_complete_list, 'variants')
        
        consensus_complete_list_noFR3, vcfs_complete_list_noFR3 = prepare_complete_consensus('-noFR3', list_spec_noFR3_bams,
                                                                                  consensus_complete,
                                                                                  folder_completevcfs,
                                                                                  rearrangement_refs)
        vcf_parsing(results_folder, folder_completevcfs, vcfs_complete_list_noFR3, 'variants')
    
    ## Consensus sequence and local alignment against the reference alleles using EMBOS water
    #water_list = mutationalStatus(consensus_complete_list, homology_folder, ref_dict)
    water_list2 = mutationalStatus(consensus_complete_list_noFR3, homology_folder, ref_dict)
    ## Parsing homology data and adding the value to dictionary
    ## We will use water_list to open each alignment file and parse EMBOSS
    ## water output to add the homology result to the final table
    #d_hom, g_hom = homology_parsing(water_list, d_hom, g_hom, 'homology', 'mutational_status')
    d_hom, g_hom = homology_parsing(water_list2, d_hom, g_hom, 'homology_noFR3', 'mutational_status_noFR3')

    prepare_consensus_sequence_annotation(d_hom, g_hom, consensus_complete_list_noFR3)

    # Productivity and CDR3
    #fasta_IGHD_list, junctiond, d_hom, g_hom = motifs(consensus_complete_list,
    #                                                  d_hom, g_hom, CDR3_fastas_folder,
    #                                                  folder_completevcfs)
    if not clonal:
        # get unique reads from bam_specific_rearrangements
        list_uniq_reads, d_hom, DH_list = bams_specific_read(list_spec_bams, uniq_bam_sequences,
                                                    folder_completevcfs, homology_folder,
                                                    ref_dict, ref_dictJ, d_hom, info_folder)

    # BLAST CDR3
    #blast_list = prepare_CDR3_blast(fasta_IGHD_list, fblast_CDR3)
    emboss_list, EW_list = prepare_CDR3_emboss(DH_list, fblast_CDR3, refD)

    ## add CDR3 blast results to homology table
    #blast_annotation(blast_list, d_hom, g_hom)
    emboss_annotation(EW_list, d_hom, g_hom)

    ## output
    # Write dictionary to out file 'homology_table.csv'
    for k in sorted(d_hom):
        if clonal:
            s = '_'.join(k.split('_')[0:2])
            v = k.split('_')[2]
            
        else:
            s = k.split('_')[0]
            v = k.split('_')[1]
        
        l = dhom_write(s, v, k, d_hom)
        f_hom.write(l + '\n')


    f_hom.close()
    if not clonal:
        for k in sorted(g_hom):
            s = k.split('_')[0]
            v = k.split('_')[1]

            l = dhom_write(s, v, k, g_hom)

            m_hom.write(l + '\n')

        m_hom.close()

    ## Open f_hom and filter it to f_hom2 using
    ## new_regions if --basal option was given
    if args.basal and not clonal:
        fone = read_file_simply(f_homology)
        for e in fone:
            if e.startswith('sample_name'):
                f_hom2.write(e + '\n')
            else:
                for l in new_regions_listV:
                    if e.startswith(l):
                        f_hom2.write(e + '\n')

        f_hom2.close()


if __name__== "__main__":
    main()

#! /usr/env/python

"""Calculate coverage for each rearrangement in the VDJ region"""

# import modules
import sys
import string
import re
import os
import argparse
import glob
import time
import inspect
import pipeline
import subprocess
## plotting
import matplotlib
import pandas as pd
import seaborn as sns
#sns.set(style="whitegrid", font_scale=2)
#sns.set_style("whitegrid")
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# import submodules
repo = os.environ['scripts_repo']
log_path = os.path.join(repo, 'NGS-tools')
sys.path.append(log_path)
from log_manager import log_setup
from dir_management import *

### FUNCTIONS ###
def extract_args():

    """"""
    # global
    global out_folder
    global wd
    global target_table
    global proc
    global ow
    global min_cov
    

def parallel(CMD, p, LOG): ## make independent function with execute

    """"""
    parallel = ('parallel --compress --joblog {} '
                '-j {} :::: {}').format(LOG, p, CMD)

    execute(parallel)
                

def execute(CMD):

    """"""
    subprocess.call([CMD], shell=True, executable='/bin/bash')

    
def create_dict_list(tab, field=2):

    """"""
    d = {}
    ##ori##
    #d = { (m.split(',')[1] + '_' + m.split(',')[2]) : ({'nreads_allele': m.split(',')[3], \
    #                                               'percent_reads_allele': m.split(',')[4],\
    #                                               'gene': m.split(',')[6], 'nreads_gene': m.split(',')[7], \
    #                                                'percent_reads_gene': m.split(',')[8], 'mut_status' : m.split(',')[10]}) for m in tab }

    d = { (m.split(',')[0] + '_' + m.split(',')[1]) : ({'nreads_allele': m.split(',')[field]}) for m in tab }
    return d


def prepare_bedgraph(bam, fcmd, out_folder, ow):

    """"""
    pipeline.owrt=ow
    overwrite = pipeline.do_overwrite(out_folder, '.bedGraph')
    if overwrite:
        log.debug('Starting bedgraph generation. Out folder %s', out_folder)
        path_bedgraph = os.path.join(out_folder, os.path.basename(bam).replace('.bam', '.bedGraph'))
        fcmd.write('bedtools genomecov -bg -ibam {} > {}\n'.format(bam, path_bedgraph))


def prepare_cov(bam, fcmd, out_folder, list_beds, ow, wd):

    """"""
    pipeline.owrt=ow
    overwrite = pipeline.do_overwrite(out_folder, '.cov')
    folder_refs = os.path.join(wd, 'rearrangement_refs')
    
    if overwrite:
        log.debug('Starting cov files generation. Out folder %s', out_folder)
        path_cov = os.path.join(out_folder, os.path.basename(bam).replace('.bam', '.cov'))
        ref = os.path.join(folder_refs,
                           '_'.join(os.path.basename(bam.replace('-sorted.bam',
                                                                 '.fa').replace('-noFR3','')).split('_')[1:]))
        bed_path = os.path.join(out_folder, os.path.basename(ref).replace('.fa', '.bed'))
        
        end = str(len(read_file_simply(ref)[1]) - 1)
        if not os.path.isfile(bed_path):
            with open(bed_path, 'w') as bed:
                bed.write(os.path.basename(ref).replace('.fa', '') + '\t65\t' + end + '\n')
                list_beds.append(os.path.basename(ref).replace('.fa', ''))
        fcmd.write('bedtools coverage -d -a {} -b {} > {}\n'.format(bed_path, bam, path_cov))
    else:
        list_beds = [ os.path.basename(e).replace('.bed', '') for e in glob.glob(os.path.join(out_folder, '*.bed'))]


    return list_beds
            
    
## main
def calculate_coverage(out_folder, wd, target_table, proc, ow, min_cov):

    """Main script for pipeline intro"""
                        
    ## Define logging
    #log_setup()
    
    ## create folders
    folder_bedgraphs = os.path.join(out_folder, 'bedgraphs')
    create_dir(folder_bedgraphs)
    folder_cov = os.path.join(out_folder, 'covs')
    create_dir(folder_cov)
    
    ## Read rearrangements table
    tab = read_file_simply(target_table)
    if "homology_resume" in target_table:
        field = 3
    else:
        field = 2
    d = create_dict_list(tab[1:], field) # without header

    ## calculate cov per sample from sorted indexed bams for the whole rearrangement
    path_bams = os.path.join(wd, 'bams_specific_rearrangements')
    # define cdate 
    pipeline.cdate=time.strftime("%d-%m-%Y_%H-%M-%S")
    path_CMD_bedgraph, path_LOG_bedgraph = pipeline.create_CMD('bedGraph', out_folder)
    path_CMD_cov, path_LOG_cov = pipeline.create_CMD('cov', out_folder)
    fcmd_bedgraph = open(path_CMD_bedgraph, 'w')
    fcmd_cov = open(path_CMD_cov, 'w')
    rearrangements_list = []
    for sample_id in d:
        sample_id = sample_id.replace('-D0', 'D-0')
        
        bam = glob.glob(path_bams + '/' + sample_id + '*.bam')
        if len(bam) > 0:
            ## calculate coverage
            prepare_bedgraph(bam[0], fcmd_bedgraph, folder_bedgraphs, ow)
            rearrangements_list = prepare_cov(bam[0], fcmd_cov,
                                              folder_cov, rearrangements_list,
                                              ow, wd)
        
    fcmd_bedgraph.close()
    fcmd_cov.close()
    parallel(path_CMD_bedgraph, proc, path_LOG_bedgraph)
    parallel(path_CMD_cov, proc, path_LOG_cov)

    ### MEAN ALL FILES
    ## FOF files
    bg_fof = glob.glob(folder_bedgraphs + '/*.bedGraph')
    
    fofbg = pipeline.list2file(bg_fof, 'bedGraphs.fof', folder_bedgraphs)
    cov_fof = glob.glob(folder_cov + '/IGH-*.cov')
    
    fofcov = pipeline.list2file(cov_fof, 'cov.fof', folder_cov)

    ## cov stats per sample
    path0 = os.path.join(out_folder, 'samplecovstats.csv')
    execute('python /srv/dev/IGH_programming/meanpersample.py {} {} > {}'.format(fofcov, min_cov, path0))
    ## python bedgraphs mean
    path1 = os.path.join(folder_bedgraphs, 'all.bedgraph')
    path2 = os.path.join(folder_cov, 'all.cov')
    execute('python /srv/dev/IGH_programming/mean_bedgraph.py {} > {}'.format(fofbg, path1))
    execute('python /srv/dev/IGH_programming/mean_covperbase_freebayes.py {} > {}'.format(fofcov, path2))

    ## PLOTS PER SAMPLE AND REARRANGEMENT
    for f in cov_fof:
        plt.figure(figsize=(25,14))
        sns.set(style="whitegrid", font_scale=2)
        d = pd.DataFrame()
        ld = read_file_simply(f)
        l1 = []
        l2 = []
        l3 = []
        for m in ld:
            name, start, end , pos, cov = m.split('\t')
            l2.append(pos)
            l3.append(cov)
    
        d['positions'] = l2
        d['coverage'] = l3
        d.coverage = d.coverage.astype(float)
        d.positions = d.positions.astype(int)
        g = sns.lineplot(x='positions', y='coverage', data=d, lw=5)
        ax1 = g.axes
        ax1.grid(False)
        ax1.axhline(int(min_cov), ls='--', c='r')
        ax1.set_title('Coverage per base {}'.format(os.path.basename(f).replace('.cov', '')))
        ax1.set(xlabel='Base', ylabel='Coverage Depth')
        plt.savefig(f.replace('.cov', '.png'))
        plt.gcf().clear()

    ### MEAN PER REARRANGEMENT
    bedgraph_per_rearrangement = os.path.join(out_folder, 'meanperrearrangement_bedgraphs')
    create_dir(bedgraph_per_rearrangement)
    cov_per_rearrangement = os.path.join(out_folder, 'meanperrearrangement_covs')
    create_dir(cov_per_rearrangement)
    
    ## create dataframe for general plot per rearrangement
    df = pd.DataFrame()
    
    for b in rearrangements_list:
    
        bg_fof = glob.glob(folder_bedgraphs + '/*' + b + '*.bedGraph')
        fofbg = pipeline.list2file(bg_fof, 'bedGraphs_' + b + '.fof', folder_bedgraphs)
        cov_fof = glob.glob(folder_cov + '/*' + b + '*.cov')
        fofcov = pipeline.list2file(cov_fof, 'cov_' + b + '.fof', folder_cov)

        path1 = os.path.join(bedgraph_per_rearrangement, 'all_' + b + '.bedgraph')
        path2 = os.path.join(cov_per_rearrangement, 'all_' + b + '.cov')

        execute('python /srv/dev/IGH_programming/mean_bedgraph.py {} > {}'.format(fofbg, path1))
        execute('python /srv/dev/IGH_programming/mean_covperbase_freebayes.py {} > {}'.format(fofcov, path2))

        # range for plot
        
        # read mean cov file
        list_r = read_file_simply(path2)
        list_names = []
        list_cov = []
        list_pos = []
        for e in list_r:
            name, pos, cov = e.split('\t')
            list_names.append(name)
            list_cov.append(cov)
            list_pos.append(pos)
        
        if df.empty:
            df['positions'] = list_pos
            df[b] = list_cov
            
        else:
            df2 = pd.DataFrame()
            df2['positions'] = list_pos
            df2[b] = list_cov
            df = pd.merge(df, df2, on='positions', how='outer')
                  
    df = df.melt('positions', var_name='cols', value_name='vals')
    
    df.positions = df.positions.astype(int)
    df.vals = df.vals.astype(float)
    df.to_csv(os.path.join(out_folder, 'covperrearrangement.csv'), sep=',')
    
    #g = sns.catplot(x='positions', y='vals', hue='cols', data=df)
    sns.set_style("whitegrid")
    g = sns.lineplot(x='positions', y='vals', data=df)
    ax = g.axes
    ax.grid(False)
    #ax.axhline(int(min_cov), ls='--', c='r')
    ax.set_ylim(ymin=0)
    ax.set_title('Mean coverage depth per rearrangement')
    ax.set(xlabel='Base', ylabel='Coverage Depth')
    plt.savefig(os.path.join(out_folder, 'mean_rearrangements.png'))
    plt.gcf().clear()

    # the same but showing all the experimental data coloured by rearrangements
    g = sns.lineplot(x='positions', y='vals', hue='cols', legend=False,
                 estimator=None, lw=1.5, data=df)
    ax = g.axes
    ax.grid(False)
    ax.axhline(int(min_cov), ls='--', c='r')
    ax.set_title('Mean coverage depth per rearrangement')
    ax.set(xlabel='Base', ylabel='Coverage Depth')
    plt.savefig(os.path.join(out_folder, 'mean_rearrangements_allcolour.png'))
    plt.gcf().clear()

    # the same but showing all the experimental data without colour
    g = sns.lineplot(x='positions', y='vals',
                 estimator=None, lw=1, data=df)
    ax = g.axes
    ax.grid(False)
    ax.axhline(int(min_cov), ls='--', c='r')
    ax.set_title('Mean coverage depth per rearrangement')
    ax.set(xlabel='Base', ylabel='Coverage Depth')
    plt.savefig(os.path.join(out_folder, 'mean_rearrangements_all.png'))
    plt.gcf().clear()

    
    
if __name__=="__main__":

    main()

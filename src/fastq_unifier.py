#! /usr/bin/python3.5

import re
import os

# This script recognizes the fastq format names to use them as the input of a pipeline.
# The considerations are that FASTQ names must be inside the universal format of FASTQ
# files and that the sample name can not contain "_" in it. In this case FASTQ names
# have to be changed (for example replace "_" by "-", though "_" is already the universal
# separator for the rest of the fields.


def merge(sample_name_read_dict, trim_dict, merged_folder, trimming, pairs_raw, d, pairs_trimmed):

    """distribution of FASTQ files per sample name and read to perform merge and trimming"""
                                                                                                                       
    merged = []
    for sample in sample_name_read_dict:
        trim_dict[sample] = {'merged_name':'', 'trimmed_name':'', 'fastqs_to_merge':''}
        trim_dict[sample]['fastqs_to_merge'] = sorted(sample_name_read_dict[sample]['name'].split(';'))
        indext = sample_name_read_dict[sample]['name'].split(';').index(trim_dict[sample]['fastqs_to_merge'][0])
        ori_fastq_path = trim_dict[sample]['fastqs_to_merge'][0]
        extension = sample_name_read_dict[sample]['extension'].split(';')[indext]
        
        trim_dict[sample]['merged_name'] = os.path.join(merged_folder, os.path.basename(ori_fastq_path))
        
        if len(os.path.basename(trim_dict[sample]['merged_name']).split('_')) <= 1:
            trim_dict[sample]['trimmed_name'] = os.path.basename(trim_dict[sample]['merged_name'].replace(('.' + extension),'_001-trimmed.fastq.gz'))
        else:
            trim_dict[sample]['trimmed_name'] = os.path.basename(trim_dict[sample]['merged_name'].replace(('.' + extension),'-trimmed.fastq.gz'))
        ## merged list
        merged.append(';'.join(trim_dict[sample]['fastqs_to_merge']))
        trimming.append(trim_dict[sample]['merged_name'] + ';' + trim_dict[sample]['trimmed_name'])

        if len(sample.split('_')) > 1:
            if '1' in sample.split('_')[1]:

                pairs_raw.append(trim_dict[sample]['merged_name'] + ';' + os.path.join(merged_folder, os.path.basename(d[ori_fastq_path]['mate'])))
                s2 = (os.path.basename(d[ori_fastq_path]['mate'])).replace(('.' + extension),'-trimmed.fastq.gz')
                pairs_trimmed.append(trim_dict[sample]['trimmed_name'] + ';' + s2)
        else:
            pairs_raw.append(trim_dict[sample]['merged_name'])
            pairs_trimmed.append(trim_dict[sample]['trimmed_name'])
            
        
    return merged, pairs_raw, pairs_trimmed, trimming, trim_dict 


def fastq_dictionary(n, merged_folder):
    
    """must provide this function with a list including all FASTQ names"""

    d = {}
    extensions = []
    pairs_raw = []
    pairs_trimmed = []
    trimming = []
    trim_dict = {}
    sample_name_dict = {}
    sample_name_read_dict = {}
    
    if '' in n:
        n.remove('')
    for fastq_path in n:
        
        # Extract basename (it will be the key of our dictionary)
        fastq = os.path.basename(fastq_path)

        # Define the empty dictionary
        d[fastq_path] = { 'path' : '', 'compressed' : '', 'name_without_extension' : '', 'extension': '',
                          'sample_name' : '', 'pair' : '', 'mate' : '', 'read' : '', 'lane' : '', 'barcode' : '',
                          'bam_name' : ''}
        
        # Extract the path of the file
        d[fastq_path]['path'] = os.path.dirname(fastq_path)

        # Manage extension and compression
        if fastq.endswith(".gz"):
            d[fastq_path]['compressed'] = 'True'
            name_without_extension = '.'.join(fastq.split('.')[:-2])
            extension = '.'.join(fastq.split('.')[-2:])
        else:
            d[fastq_path]['compressed'] = 'False'
            name_without_extension = '.'.join(fastq.split('.')[:-1])
            extension = '.'.join(fastq.split('.')[-1:]) 

        # Keep name without extension
        d[fastq_path]['name_without_extension'] = name_without_extension
        d[fastq_path]['extension'] = extension
        extensions.append(extension)
        name_list = name_without_extension.split('_')
        # We expect the first element of our list to be the sample name
        sample_name = name_list[0]
        d[fastq_path]['sample_name'] = sample_name
        d[fastq_path]['bam_name'] =  sample_name + '.bam'
        d[fastq_path]['sortedbam_name'] = sample_name + '-sorted.bam'
        name_list_without_samplename = name_list[1:]      

        # trimmed_name
        if len(name_list_without_samplename) > 0:
            d[fastq_path]['trimmed_name'] = sample_name + '_' + '_'.join(name_list_without_samplename) + '-trimmed.' + extension
        else:
            d[fastq_path]['trimmed_name'] = sample_name + '_001-trimmed.' + extension
            
        # detect read (R1, R2)
        paired = True
        if 'R1' in name_list_without_samplename:
            inde = name_list_without_samplename.index('R1')
            pair = 'forward'
            mate = (sample_name + '_' + '_'.join(name_list_without_samplename)).replace('_R1', '_R2') + '.' + extension
                     
        elif 'R2' in name_list_without_samplename:
            inde = name_list_without_samplename.index('R2')
            pair = 'reverse'
            mate = (sample_name + '_' + '_'.join(name_list_without_samplename)).replace('_R2', '_R1') + '.' + extension
            
            
        elif '1' in name_list_without_samplename:
            inde = name_list_without_samplename.index('1')
            pair = 'forward'
            mate = (sample_name + '_' + '_'.join(name_list_without_samplename)).replace('_1', '_2') + '.' + extension
                        
        elif '2' in name_list_without_samplename:
            inde = name_list_without_samplename.index('2')
            pair = 'reverse'
            mate = (sample_name + '_' + '_'.join(name_list_without_samplename)).replace('_2', '_1') + '.' + extension
                    
        else:
            paired = False
            
        if paired:
            d[fastq_path]['pair'] = pair
            d[fastq_path]['mate'] = os.path.join(d[fastq_path]['path'], mate)
            read = name_list_without_samplename[inde]
            d[fastq_path]['read'] = read
                              
        if d[fastq_path]['sample_name'] not in sample_name_dict:
            sample_name_dict[d[fastq_path]['sample_name']] = fastq_path
        else:
            sample_name_dict[d[fastq_path]['sample_name']] = sample_name_dict[d[fastq_path]['sample_name']] + ';' + fastq_path
        
        # create dictionary per sample name and read
        if paired:
            identifier = sample_name + '_' + read
        else:
            identifier = sample_name
            
        if not identifier in sample_name_read_dict:
            sample_name_read_dict[identifier] = {'name' : '', 'extension' : ''}
            sample_name_read_dict[identifier]['name'] = fastq_path
            sample_name_read_dict[identifier]['extension'] = d[fastq_path]['extension']
        else:
            sample_name_read_dict[identifier]['name'] = sample_name_read_dict[identifier]['name'] + ';' + fastq_path
            sample_name_read_dict[identifier]['extension'] = sample_name_read_dict[identifier]['extension'] + ';' + d[fastq_path]['extension']
    
        # detect lane and barcode (6 nucleotides)
        for element in name_list_without_samplename:
            if element.startswith('L'):
                lane = element
                d[fastq_path]['lane'] = lane
            elif re.match('(?=[TCGA]){6}', element):
                barcode = element
                d[fastq_path]['barcode'] = barcode

    merged, pairs_raw, pairs_trimmed, trimming, trim_dict = merge(sample_name_read_dict, trim_dict,
                                                                  merged_folder, trimming, pairs_raw, d, pairs_trimmed)
                                                            
    return (d, extensions, sample_name_dict, sample_name_read_dict, merged,
            pairs_raw, set(pairs_trimmed), trimming, trim_dict)

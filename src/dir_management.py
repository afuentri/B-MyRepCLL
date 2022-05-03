# ! /usr/env/python3.5                      #
# Azahara Maria Fuentes Trillo              #
# Unidad de Genomica y Diagnostico Genetico #
# INCLIVA                                   #
# 2018-04-11                                #
#############################################

# Import modules
import os
import string
import gzip
import logging as log
import subprocess
import re
import sys

# logging

def check_directory(dire):

    """Check if a given path is a directory and if not it creates it"""
    if not os.path.isdir(dire):
        try:
            os.makedirs(dire)
        except OSError as e:
            sys.exit("Error creating {} directory. Check path and permissions".format(dire))
        except:
            sys.exit("Unexpected error:", sys.exc_info()[0])
            raise
            

def check_empty(dire):
    """Check if a device is empty or not mounted"""
    if os.path.isdir(dire):
        if os.listdir(dire) == []:
            log.critical('%s is empty, exiting script', dire)
        
        else:
            log.debug('%s is correctly mounted', dire)
    else:
        log.warning('%s is not a directory', dire)
        
def check_path_exists(path):

    """Check if the given path exists"""
    if not os.path.exists(path):
        log.error("The path %s does not exist", path)
        

def open_file(input, mode='read'):

   """Open file if exists"""
   while True:
       try:
           if mode == 'write':
               f_hand = open(input, 'w')
           elif mode == 'read':
               f_hand = open(input, 'r')
           elif mode == 'add':
               f_hand = open(input, 'a')
       except:
           raise IOError(input + ' is not a file. Try again')
           continue
       else:
           log.info('Reading %s successfully', input)
           break


   return f_hand


def read_file_simply_gzip(input):

    """Read a file compressed with gzip and return a list with each line as element"""
    while True:
        try:
            f_hand = gzip.open(input, 'r')
        except:
            raise IOError(input + ' is not a file. Try again')
            continue
        else:
            log.info('Reading %s successfully', input)
            with gzip.open(input, 'rt', encoding='utf-8') as f_hand:
                list = f_hand.read().splitlines()

            f_hand.close()
            break

    return list

                
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
            log.info('Reading %s successfully', input)
            with open(input, encoding='utf-8') as f_hand:
                list = f_hand.read().splitlines()

            f_hand.close()
            break

    return list

                            
def create_dir(directory):

    """Create a directory"""
    if not os.path.exists(directory):
        os.makedirs(directory)


# ! /usr/env/python3.5                      #
# Azahara Maria Fuentes Trillo              #
# Unidad de Genomica y Diagnostico Genetico #
# INCLIVA                                   #
# 2018-04-11                                #
#############################################
import sys
import logging as log
from logging import handlers

# LOGGING LEVELS
# DEBUG Detailed information, typically of interest only when diagnosing problems.
# INFO Confirmation that things are working as expected.
# WARNING An indication that something unexpected happened, or indicative of some problem in the near future (e.g. 'disk space low'). The software is still working as expected.
# ERROR Due to a more serious problem, the software has not been able to perform some function.
# CRITICAL A serious error, indicating that the program itself may be unable to continue running.

def log_setup(args):

    """Define logging mode"""
    # Logging and verbose
    if args.verbose:
        log.basicConfig(format='%(levelname)s: %(message)s', level=log.DEBUG)
    else:
        log.basicConfig(format='%(levelname)s: %(message)s', level=log.INFO)

    #logger = log.getLogger(__name__)

    # Create a logging handler for loglevel warning
    handler = log.FileHandler('warnings.log')
    handler.setLevel(log.INFO)
    ## file handler logging format
    formatter = log.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    ## add the handlers to the logger
    #log.addHandler(handler)
    log.info('Starting warning handling to file warnings.log')


    

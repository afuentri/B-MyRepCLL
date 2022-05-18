import sys
import glob
import os
import subprocess
import string
import time

proj_name = sys.argv[1]
mincov = sys.argv[2]
folder_fastqs = sys.argv[3]
primers = sys.argv[4]

"""
Simple Python script with the only purpose of launching 
IGH pipeline.
"""

def execute(CMD):
    
    """"""
    subprocess.call([CMD], shell=True, executable='/bin/bash')

## variables
repo = os.path.abspath(__file__ + "/../../")
cdate = time.strftime("%Y%m%d")
nfastqs = len(glob.glob(str(folder_fastqs) +
                        '/*fastq.gz'))

phrase = 'Detected {} FASTQ files for {} samples. Running on {} processes.'.format(nfastqs, nfastqs/2, '20')
CMD2 = ('time python {3}/B-MyRepCLL/src/pipeline.py '
        '--pipeline -f {0} '
        '-o {1}_{2} -v -p15 --basal '
        '--primers {4} '
        '--cdr3s > {1}_{2}.log').format(folder_fastqs, cdate, proj_name, repo, primers)

execute(CMD2)

## QC
QCMD = ('time bash {3}/QC/main-parallel.sh -p {0}_{1} '
        '-b {2} -t 10').format(cdate, proj_name, os.path.dirname(primers), repo)

execute(QCMD)

## flagstat
FCMD = ('time bash {2}/QC/flagstat.sh -b '
        '{0}_{1}/bamsV').format(cdate, proj_name, repo)

execute(FCMD)

## filter and merge artifactual rearrangements and generate final summary files
CMD1 = ('time python3 {2}/B-MyRepCLL/src/filter_equalrearrangements.py '
        '{0}_{1}/results/homology_table.csv '
        '{0}_{1}/results/consensus_complete').format(cdate, proj_name, repo)

execute(CMD1)

CMD2 = ('time python3 {2}/B-MyRepCLL/src/filter_equalrearrangements2.py '
        '{0}_{1}/results/homology_tablesimpalleles.csv '
        '{0}_{1}/results/consensus_complete').format(cdate, proj_name, repo)

execute(CMD2)

CMD3 = ('time python3 {3}/B-MyRepCLL/src/onlyclonality.py '
        '{1}_{0}/results/homology_tablesimpallelessimpalleles.csv '
        '{0} {1}_{0} {1}_{0}/QC/fastq_stats.xlsx '
        '{1}_{0}/bamsV/flagstat/resume.csv {2}').format(proj_name, cdate, mincov, repo)

execute(CMD3)

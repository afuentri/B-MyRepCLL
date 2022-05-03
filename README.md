# B-MyRepCLL
pipeline for a consensus-oriented analysis of B-cell clones

The strategy is mapping reads separately against the different IMGT gene segments references, and following a clone-centered determination, achieved with the obtention of a consensus sequence. B-cell rearrangements are defined after IGHV-IGHJ correspondence determination and a specific procedure has been designed to cope with unspecific mapping and gene-call primer biases, and the calculation of the clonal fraction  the unique profile of each patient. 

The first module, src/pipeline.py, annotates VDJ calls and mutational status for all the IGHV alleles found per patient. The second module, src/onlyclonality.py, generates filtering steps for the minimization of artifacts and outputs homology_resume*.xlsx with the final results, with a calculation of the clonal and subclonal fraction on each sample.

![alt text](pipeline.png)

For installation, first clone the repository.

You can install all the dependencies directly using conda:
```console
conda env create -f environment.yml
```

if not, install the following requirements manually:
### Requirements
>bwa 0.7.15 or above  
>bamtools 2.4 or above  
>bcftools 1.7 or above  
>bedtools 2.26 or above  
>bbduk(bbtools), repair(bbtools) BBMap version 38 or above  
>emboss water 6.6.0 or above  
>samtools 1.7 or above  
>freebayes 1.1.0 or above  
>seqtk 1.2 or above  
>Python  
>R  

#### PIP Install Requirements
```console
## Python 2
pip install -r requirements.txt 
## or Python 3
pip3 install -r requirements.txt
```

#### Recommended arguments (default workflow validated with CLL data)
```console
python B-MyRepCLL/src/pipeline.py --pipeline -f $fastqfilesFolder -o $outputDir -v -p$nproc --basal --primers $FASTAfileVHprimers --cdr3s > log.log
```
#### Automatic execution of the pipeline with default validated parameters, final summary files and quality control
```console
python B-MyRepCLL/launch-default.py $fastqfilesFolder $coverage_threshold $outputDir
```
This mode has requirements of other repositories:
https://github.com/afuentri/QC

Publication in process

bwa mem -t4 -M -R "@RG\tID:IGH-MIX3-clonal-control\tPL:ILLUMINA\tSM:IGH-MIX3-clonal-control\tDS:ref=/nfs/hpcugdga/IGHS-motivos_20171205/references/IMGT_J-nostars-simple.fa\tCN:UGDG\tDT:15-10-2019_20-18-16\tPU:IGH-MIX3-clonal-control" /nfs/hpcugdga/IGHS-motivos_20171205/references/IMGT_J-nostars-simple.fa /media/scratch/aft_projects_tmp/python3_myenv/test/trimmed/IGH-MIX3-clonal-control_S58_L001_R1_001-trimmed.fastq.gz /media/scratch/aft_projects_tmp/python3_myenv/test/trimmed/IGH-MIX3-clonal-control_S58_L001_R2_001-trimmed.fastq.gz | samtools view -b - | samtools sort -o /media/scratch/aft_projects_tmp/python3_myenv/test/bamsJ/IGH-MIX3-clonal-control-sorted.bam -

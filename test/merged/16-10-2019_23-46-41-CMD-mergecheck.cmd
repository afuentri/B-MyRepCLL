n1=$(zcat /media/scratch/aft_projects_tmp/pigh/test/fastq_test/IGH-MIX3-clonal-control_S58_L001_R1_001.fastq.gz | wc -l); n2=$(zcat /media/scratch/aft_projects_tmp/pigh/test/merged/IGH-MIX3-clonal-control_S58_L001_R1_001.fastq.gz | wc -l); if [ $n1 -eq $n2 ]; then echo "Merge was correct for sample /media/scratch/aft_projects_tmp/pigh/test/merged/IGH-MIX3-clonal-control_S58_L001_R1_001.fastq.gz"; else echo "/media/scratch/aft_projects_tmp/pigh/test/merged/IGH-MIX3-clonal-control_S58_L001_R1_001.fastq.gz" >> /media/scratch/aft_projects_tmp/pigh/test/merged/16-10-2019_23-46-41-wrongmerge.txt; fi
n1=$(zcat /media/scratch/aft_projects_tmp/pigh/test/fastq_test/IGH-MIX3-clonal-control_S58_L001_R2_001.fastq.gz | wc -l); n2=$(zcat /media/scratch/aft_projects_tmp/pigh/test/merged/IGH-MIX3-clonal-control_S58_L001_R2_001.fastq.gz | wc -l); if [ $n1 -eq $n2 ]; then echo "Merge was correct for sample /media/scratch/aft_projects_tmp/pigh/test/merged/IGH-MIX3-clonal-control_S58_L001_R2_001.fastq.gz"; else echo "/media/scratch/aft_projects_tmp/pigh/test/merged/IGH-MIX3-clonal-control_S58_L001_R2_001.fastq.gz" >> /media/scratch/aft_projects_tmp/pigh/test/merged/16-10-2019_23-46-41-wrongmerge.txt; fi

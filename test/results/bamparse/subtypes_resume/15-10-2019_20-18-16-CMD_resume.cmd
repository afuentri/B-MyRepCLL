cat /media/scratch/aft_projects_tmp/python3_myenv/test/results/bamparse/bamparsing_out/info_bams_IGH-MIX3-clonal-control.txt | grep "region" | cut -f2,4 | sort | uniq -c | sort -rn > /media/scratch/aft_projects_tmp/python3_myenv/test/results/bamparse/subtypes_resume/subtypes_resume_IGH-MIX3-clonal-control.txt

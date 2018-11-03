#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o ../out/logs                        #-- output directory (fill in)
#$ -e ../out/logs                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=128G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=25G,scratch=1G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=11:58:00                #-- runtime limit (see above; this requests 24 hours)
#$ -M wconnell93@gmail.com
#$ -m bea

date
hostname

qstat -j $JOB_ID

echo $input
Rscript monocle_production.R $input

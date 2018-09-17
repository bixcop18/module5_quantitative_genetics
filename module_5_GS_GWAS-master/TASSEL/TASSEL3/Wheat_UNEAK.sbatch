#!/usr/bin/env bash


###########################################
## Shell script to run GBS in TASSEL 3   ##
##             UNEAK Pipeline            ##
###########################################

##These commands seem to cause problems but could be used as follows:
##SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END 
##SBATCH --mem-per-cpu=30G   # Memory per core, use --mem= for memory per node
#SBATCH -p batch
#SBATCH -n 1
#SBATCH --time=0-06:00:00   # Use the form DD-HH:MM:SS
#SBATCH -J wheat_example  #Add a name to the job

WORKDIR=/var/scratch/$USER/$SLURM_JOBID
mkdir -p $WORKDIR

echo "Using $WORKDIR on $SLURMD_NODENAME"
echo

# change to working directory on compute node
cd $WORKDIR

#load modules
module load java/1.8  #load java

#Give analysis a Name
NAME=wheat_example


TASSEL3=/home/jcrain/TASSEL_Software/tassel3-standalone/run_pipeline.pl


#Step 1 create working directories
${TASSEL3} -Xms30G -Xmx30G -fork1 -UCreatWorkingDirPlugin  -w ${WORKDIR} -endPlugin -runfork1

#Step 2
#copy fastq file and key into respective directories
cp /home/jcrain/examples/wheat/fastq/*.txt.gz $WORKDIR/Illumina/
cp /home/jcrain/examples/wheat/keys/SynOpDH_UNEAK_Key_Multiple_Fastq.txt $WORKDIR/key/


#Step 3 Make tag count file
${TASSEL3} -Xms30G -Xmx30G -fork1 -UFastqToTagCountPlugin -w ${WORKDIR} -e PstI-MspI -s 500000000 -c 1 -endPlugin -runfork1

#Step 4 Make master tag count file 
${TASSEL3} -Xms30G -Xmx30G -fork1 -UMergeTaxaTagCountPlugin -w ${WORKDIR} -c 5  -m 1000000000 -x 10000000 -endPlugin -runfork1

#Step 5
${TASSEL3} -Xms30G -Xmx30G -fork1 -UTagCountToTagPairPlugin -w ${WORKDIR} -e 0.03 -endPlugin -runfork1

########Generate Tags by Taxa TBT file###############     
${TASSEL3} -Xms30G -Xmx30G -fork1 -UTagPairToTBTPlugin -w ${WORKDIR} -endPlugin -runfork1

######Generate map info file #####################
${TASSEL3} -Xms30G -Xmx30G -fork1 -UTBTToMapInfoPlugin -w ${WORKDIR} -endPlugin -runfork1
     
##########Generate HapMap file ###############
${TASSEL3} -Xms30G -Xmx30G -fork1 -UMapInfoToHapMapPlugin -w ${WORKDIR} -mnMAF 0.05 -mxMAF 0.5     -mnC 0  -mxC 1 -endPlugin -runfork1

#Make sure to manually cleanup /var/scratch/${user}.  Delete all binary files and transfer hapmap file to home directory.



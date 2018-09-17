#!/usr/bin/env bash


###########################################
## Shell script to run GBSv2 in TASSEL 5 ##
## Adapted from Narinder Singh           ##
###########################################

##These commands seem to cause problems but could be used as follows:
##SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END 
#SBATCH --mem-per-cpu=20G   # Memory per core, use --mem= for memory per node
#SBATCH -p batch
#SBATCH -n 1
#SBATCH --time=0-06:00:00   # Use the form DD-HH:MM:SS
#SBATCH --J rice_example  #Add a name to the job

WORKDIR=/var/scratch/$USER/$SLURM_JOBID
mkdir -p $WORKDIR

echo "Using $WORKDIR on $SLURMD_NODENAME"
echo

# change to working directory on compute node
cd $WORKDIR

#load modules
module load bowtie2/2.3.4.1  #load bowtie

module load java/1.8  #load java

#Give analysis a Name
NAME=rice_example


#put path to keyfile
KEYFILE=/home/jcrain/examples/rice/keys/D0DTNACXX_1_key_file_Tassel_V.tsv 

#put directory to sequence file
SEQUENCE=/home/jcrain/examples/rice/fastq

######Run one time if needed##########
######If ran once can simply point to reference####
#########Build reference genome
bowtie2-build /home/jcrain/examples/rice/reference/Osative_204_Reference.fa ${WORKDIR}/Osative_204_Reference

#put path to reference
REFERENCE=${WORKDIR}/Osative_204_Reference

#put path to tassel5 pipeline 
TASSEL=/home/jcrain/TASSEL_Software/tassel-5-standalone/run_pipeline.pl

 
## NO NEED TO CHANGE ANYTHING FROM HERE ON ##
 
#Step 1
#have to request memory but aslo specify it calling the plugin

## GBSSeqToTagDBPlugin  - RUN Tags to DB  require min quality score, 50 base pair tags, and up to 250M kmers in the database
$TASSEL -Xms20G -Xmx20G -fork1 -GBSSeqToTagDBPlugin -e ApeKI \
    -i $SEQUENCE \
    -db ${NAME}.db \
    -k ${KEYFILE} \
    -kmerLength 64 -minKmerL 50 -mnQS 20 -mxKmerNum 250000000 \
    -endPlugin -runfork1 >> ${NAME}_pipeline.out

#Step 2
## TagExportToFastqPlugin  - export Tags to align to reference
$TASSEL -fork1 -TagExportToFastqPlugin \
    -db ${NAME}.db \
    -o ${NAME}_tagsForAlign.fa.gz -c 10 \
    -endPlugin -runfork1 >> ${NAME}_pipeline.out
    
#Step 3
## RUN BOWTIE #-S is write to SAM file -U is unparied reads to be aligned -x is aligned files
bowtie2 --end-to-end  \
    -x ${REFERENCE} \
    -U ${NAME}_tagsForAlign.fa.gz \
    -S ${NAME}.sam >> ${NAME}_pipeline.out
    
#Step 4  
## SAMToGBSdbPlugin - SAM to DB, update database with alignment information
$TASSEL -Xms20G -Xmx20G -fork1 -SAMToGBSdbPlugin \
    -i ${NAME}.sam \
    -db ${NAME}.db \
    -aProp 0.0 -aLen 0 \
    -endPlugin -runfork1 >> ${NAME}_pipeline.out

#Step 5
## DiscoverySNPCaller
$TASSEL -Xms20G -Xmx20G -fork1 -DiscoverySNPCallerPluginV2 \
    -db ${NAME}.db \
    -mnLCov 0.1 -mnMAF 0.01 -deleteOldData true \
     -endPlugin -runfork1 >> ${NAME}_pipeline.out
  
#Step 6  
## SNPQualityProfilerPlugin - RUN QUALITY PROFILER
$TASSEL -Xms20G -Xmx20G -fork1 -SNPQualityProfilerPlugin \
    -db ${NAME}.db \
    -statFile ${NAME}_SNPqual_stats.txt \
    -endPlugin -runfork1 >> ${NAME}_pipeline.out
  
#Step 7    
## UpdateSNPPositionQualityPlugin - UPDATE DATABASE WITH QUALITY SCORE fast < 30 minutes 15GB
$TASSEL -Xms20G -Xmx20G -fork1 -UpdateSNPPositionQualityPlugin \
    -db ${NAME}.db \
    -qsFile ${NAME}_SNPqual_stats.txt \
    -endPlugin -runfork1 >> ${NAME}_pipeline.out


#Ends SNP discovery with database
#Use Production SNP caller to get SNPs and filter
#Step 8    
## ProductionSNPCallerPluginV2 - RUN PRODUCTION PIPELINE - output .vcf
$TASSEL -Xms20G -Xmx20G -fork1 -ProductionSNPCallerPluginV2 \
    -db ${NAME}.db \
    -i ${SEQUENCE} \
    -k ${KEYFILE} \
    -o ${NAME}.vcf \
    -e ApeKI -kmerLength 64 \
    -endPlugin -runfork1 >>  ${NAME}_pipeline.out 
    
#Step 9 ## Convert to Hapmap format
$TASSEL -Xms20G -Xmx20G -fork1 -vcf ${NAME}.vcf \
	-export ${NAME} -exportType Hapmap >>  ${NAME}_pipeline.out 
	
	
##Step 10 #clean up
#make sure to remove unneeded files other files and compress the vcf file	
	

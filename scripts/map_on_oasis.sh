#! /bin/bash
################################################################################

### SUMMARY OF FUNCTIONS ###
# One stop destination for all your mapping needs. This script, which is meant
# to be executed on the TSCC, will create qsub scripts (and execute them) to
# perform the following:
# 1. Copy raw data (fastq.gz) from the Glassome to the TSCC (without using qsub)
#    into an automatically generated directory on the Oasis file system
# 2. Decompress the files if necessary
# 3. Map the fastq files to a reference genome, producing a *.sam file
#    * using STAR or Bowtie2
# 4. Creates tag directories
# 5. Calculates the PBC coefficient - which involves removing all non-uniquely 
#    mapped reads, duplicate reads, and calculating the pileups
#    * read more about the PBC coefficient here: 
#      http://genome.ucsc.edu/ENCODE/qualityMetrics.html
# 6. Moves all files from the TSCC back to Glassome
# 7. Cleans up if neccessary
# 8. Produces a summary of executed jobs - namely which ones failed
###


### OPTIONS AND ARGUMENTS ###
# -l map files that are on the TSCC already
# -t generate qsub scripts but do not execute them
### 

### set default options ###

map_local_files=false
testing=false
glassome_path='/projects/ps-glasslab-data/'
mappingScripts_path='/projects/ps-glasslab-bioinformatics/glassutils/mapping_scripts/'
bowtie_index_path='/projects/ps-glasslab-bioinformatics/software/bowtie2/indexes/'
star_path='/projects/ps-glasslab-bioinformatics/software/STAR/'
bowtie_path='/projects/glass-group/bioinformatics/bowtie2'

# check number of arguments
if [ $# -lt 4 ] 
then
    echo "Usage: "
    echo "map_on_oasis.sh <experiment type (chip|rna)> <genome> \
<email> <input file directory> [optional arguments]"
    echo "Options:
-l    map files already on tscc or are already copied over
-t    generate qsub scripts but do not execute them"
    exit 1
fi

### parse the input ###

OPTIND=5
while getopts "lt" option ; do # set $o to the next passed option
    case "$option" in  
    l)  
       map_local_files=true 
    ;;  
    t)  
        testing=true
    ;;  
    esac
done

experimentType=$1
genome=$2
email=$3
inputDirectory=$4
###

echo "Beginning processing for $experimentType exeriments."
echo "Data contained in $inputDirectory will be mapped to the $genome genome"
echo "Email notifications will be sent to $email"

### check arguments ###

#if [ $fileSource == "glassome" ]
if ! $map_local_files
then
    inputDirectory=$(readlink -fm ${glassome_path}/${inputDirectory/data//})
else
    if [[ $inputDirectory == "/oasis/tscc/scratch/"* ]]
    then
        inputDirectory=$(readlink -fm $inputDirectory)
    else
        inputDirectory=$(readlink -fm ${glassome_path}/${inputDirectory/data//})
    fi
fi

# check that experiment type is rna (for RNA-seq) or chip (for ChIP-seq and etc)
if [ ! $experimentType == "rna" ] && [ ! $experimentType == "chip" ]
then
    echo "Error! valid choices for experiment type are 'chip' or rna' only"
    exit 1
fi

# check that the input directory exists on specified fileSource machine
if $map_local_files
then
    if [ ! -d $inputDirectory ]
    then
        echo "Error! $inputDirectory cannot be found on the TSCC - try removing\
 the -l option" 
        exit 1
    fi    
else
    if [ ! -d $inputDirectory ]
    then
        echo "Error! $inputDirectory cannot be found on the Glassome - try \
using the -l option" 
        exit 1
    fi    
fi

###

# create directory where data will be stored on glassome
glassomeOutputDirectory="/projects/ps-glasslab-data/scratch/$USER/${inputDirectory##*/}"
if [ ! -d $glassomeOutputDirectory ]
then
    mkdir -p $glassomeOutputDirectory
else
    read -p "This script will copy output files to $glassomeOutputDirectory, \
which already exists! Would you like to delete it. \
Enter y for yes and n for no [yn]?" -n 1 -r 
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        rm -rf $glassomeOutputDirectory
        mkdir -p $glassomeOutputDirectory
    fi
fi

#### copy files to oasis ###

if ! [[ $inputDirectory == "/oasis/tscc/scratch/"* ]]
then
    outputDirectory="/oasis/tscc/scratch/$USER/${inputDirectory##*/}"
    if [ -d $outputDirectory ]
    then
        read -p "$outputDirectory already exists on tscc! \
Would you like to delete it and recopy files?\
 Enter y for yes and n for no [yn]" -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]
        then
            echo "removing $outputDirectory"
            rm -rf $outputDirectory
            echo "Copying files from $inputDirectory to $outputDirectory"
            scp -r $inputDirectory $outputDirectory
        else
            echo "Files in $inputDirectory won't be copied to $outputDirectory"
        fi
    else
        echo "Copying files from $inputDirectory to $outputDirectory"
        scp -r $inputDirectory $outputDirectory
    fi
else
    outputDirectory=$inputDirectory
fi


###

### decompress fastq.gz files

echo "Decompressing raw data (fastq.gz files)"

# find fastq.gz files
compressedDirs=()
compressedPaths=( $(find $outputDirectory -path *fastq.gz -type f) )
for f in ${compressedPaths[*]}
do
    # remove file name to get sample directory
    compressedDir=${f%/*gz}
    # append sample directory to list
    compressedDirs[${#compressedDirs[*]}]=$compressedDir
done

# filter out duplicated directories
sampleDirs=$(echo "${compressedDirs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')

for sample_dir in ${sampleDirs[*]}
do
    bname=`basename $sample_dir`
    echo "Decompressing $bname"

    # If there is exactly one fastq file, use that; otherwise...
    if ls $sample_dir/*.fastq &> /dev/null; then
        if [ `ls -l $sample_dir/*.fastq | wc -l` -ne 1 ]; then
            cat $sample_dir/*.fastq > $sample_dir/$bname.fastq_joined
            rm $sample_dir/*.fastq
            mv $sample_dir/$bname.fastq_joined $sample_dir/$bname.fastq
        fi  
        # else, only one .fastq; will be used.
    else
        # If there are any .sra files, dump to .fastq
        if ls $sample_dir/*.sra &> /dev/null; then 
            for sra in $sample_dir/*.sra
                do  
                    current_dir=`pwd`
                    # CD in so that fastq-dump works correctly
                    cd $sample_dir
                    fastq-dump $sra
                    rm $sra
                    cd $current_dir
                done
            # Then compile all the .fastq
            if [ `ls -l $sample_dir/*.fastq | wc -l` -ne 1 ]; then
                cat $sample_dir/*.fastq > $sample_dir/$bname.fastq_joined
                rm $sample_dir/*.fastq
                mv $sample_dir/$bname.fastq_joined $sample_dir/$bname.fastq
            else
                # Rename singular dumped sra file
                mv $sample_dir/*.fastq $sample_dir/$bname.fastq
            fi  
        fi  
        # Make single file, unzipping simultaneously if they are zipped
        if ls $sample_dir/*.gz &> /dev/null; then 
            zcat $sample_dir/*.gz > $sample_dir/$bname.fastq
        fi  
    fi  
done

# create output directories

# make directory for tag directories on Glassome
if [ ! -d $glassomeOutputDirectory/tag_directories ]
then
    mkdir -p $glassomeOutputDirectory/tag_directories
fi

# make directory for log files on Glassome
if [ ! -d $glassomeOutputDirectory/log_files ]
then
    mkdir -p $glassomeOutputDirectory/log_files
fi

if [ ! -d $outputDirectory/qsub_scripts ]
then
    mkdir $outputDirectory/qsub_scripts
else
    # delete existing scripts
    rm $outputDirectory/qsub_scripts/*
fi

# make directory for sam files
if [ ! -d $outputDirectory/sam_files ]
then
    mkdir $outputDirectory/sam_files
fi

# make directory for tag directories on tscc
if [ ! -d $outputDirectory/tag_directories ]
then
    mkdir $outputDirectory/tag_directories
fi

# make directory for log files on tscc
if [ ! -d $outputDirectory/log_files ]
then
    mkdir $outputDirectory/log_files
fi

# make scratch directory for pbc calculation
if [ ! -d $outputDirectory/pbc ]
then
    mkdir $outputDirectory/pbc
fi

### generate qsub scripts ###

# generate a UUID for this set of jobs
uuid=$[ 1 + $[ RANDOM % 10000 ]] # generate a random number between 0 and 10000
uuid=${USER}_${uuid}

# find directory where script is located
codebase=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) 

# generate script to map fastq to genome

# find all fastq files
for sampleDir in ${sampleDirs[*]}
    do
    fastqFile=$(readlink -fm $sampleDir/*fastq)
    currentDirectory=${fastqFile%/*fastq}
    sampleName=${fastqFile%.fastq}
    sampleName=${sampleName##/*/}
    samName=""
    logName=""

    # map file
    if [ $experimentType == "chip" ]
    then
        samName="${sampleName}.${genome}.bowtie2.sam" # change extension to sam
        logName="${sampleName}.${genome}.bowtie2.log" # remove path preceding file name

        # execute bowtie
        command="$bowtie_path/bowtie2 \
-p 8 \
-x $bowtie_index_path/$genome \
$fastqFile \
> $outputDirectory/sam_files/$samName \
2> $outputDirectory/log_files/$logName \n"

        # create tag directory
        command+="makeTagDirectory \
$outputDirectory/tag_directories/$sampleName \
-genome $genome \
-checkGC $outputDirectory/sam_files/$samName \
-format sam\n"

    elif [ $experimentType == "rna" ]
    then
        samName="${sampleName}.${genome}.star.sam" # change extension to sam
        logName="${sampleName}.${genome}.star.log" # remove path preceding file name
        # execute star
        command="$star_path/STAR \
--genomeDir $star_path/genomes/$genome \
--readFilesIn $fastqFile \
--outFileNamePrefix $currentDirectory/ \
--runThreadN 4\n"
        # rename aligned file
        command+="mv $currentDirectory/Aligned.out.sam \
$outputDirectory/sam_files/$samName\n"
        # rename log file
        command+="mv $currentDirectory/Log.final.out \
$outputDirectory/log_files/$logName\n"
        # create tag directory
        command+="makeTagDirectory \
$outputDirectory/tag_directories/$sampleName \
-genome $genome \
-checkGC $outputDirectory/sam_files/$samName \
-format sam -flip\n"

    else
        echo "Error! valid choices for experiment type are 'chip' or rna' only"
        exit 1
    fi

    # calculate PBC coefficient

    uniqueFile=$outputDirectory/pbc/${sampleName}.unique.bam
    sortedFile=$outputDirectory/pbc/${sampleName}.sorted
    pileupFile=$outputDirectory/pbc/${sampleName}.pileup

    command+="samtools view -Sbq 1 $outputDirectory/sam_files/$samName > \
$uniqueFile\n"
    command+="samtools sort $uniqueFile $sortedFile\n"
    command+="samtools mpileup ${sortedFile}.bam > $pileupFile\n"
    command+="PBC=\$(awk 'BEGIN {N1=0;ND=0} {if(\$4==1){N1+=1} ND+=1} END{print N1/ND}' ${pileupFile})\n"
    command+="echo -e \"PBC    \$PBC\" >>$outputDirectory/log_files/$logName\n" #"

    # copy log file to tag directory
    command+="cp $outputDirectory/log_files/$logName \
$outputDirectory/tag_directories/$sampleName\n"
    
    # copy files to Glassome scratch directory
    # copy log file
    command+="cp $outputDirectory/log_files/$logName \
$glassomeOutputDirectory/log_files\n"
    command+="cp -r $outputDirectory/tag_directories/$sampleName \
$glassomeOutputDirectory/tag_directories\n"
    


    # create qsub script
    echo -e "#!/bin/bash
#PBS -q hotel
#PBS -N ${sampleName}_${experimentType}_${genome}
#PBS -l nodes=1:ppn=8
#PBS -l walltime=4:00:00
#PBS -o $outputDirectory/qsub_scripts/${sampleName}_torque_output.txt
#PBS -e $outputDirectory/qsub_scripts/${sampleName}_torque_error.txt
#PBS -V
#PBS -M $email
#PBS -m abe
#PBS -A glass-group
$command" > $outputDirectory/qsub_scripts/${sampleName}.torque.sh

    # submit script
    if ! $testing
        then
        echo "Submitting job for $sampleName"
        qsub $outputDirectory/qsub_scripts/${sampleName}.torque.sh
    fi
done

### LICENSE STATEMENT ###
# Copyright (c) 2015, Jenhan Tao 
# All rights reserved. 
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met: 
#
# * Redistributions of source code must retain the above copyright notice, 
#   this list of conditions and the following disclaimer. 
# * Redistributions in binary form must reproduce the above copyright 
#   notice, this list of conditions and the following disclaimer in the 
#   documentation and/or other materials provided with the distribution. 
# * Neither the name of UC San Diego nor the names of its contributors may 
#   be used to endorse or promote products derived from this software 
#   without specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE. 
###
###############################################################################

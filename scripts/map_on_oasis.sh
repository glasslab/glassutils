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
map_only=false
no_emails=false
paired=false
copy_sam=false
glassome_path='/projects/ps-glasslab-data/'
scratch_path='/projects/ps-glasslab-scratch/'
mappingScripts_path='/projects/ps-glasslab-bioinformatics/glassutils/mapping_scripts/'
scripts_path='/projects/ps-glasslab-bioinformatics/glassutils/scripts/'
bowtie_index_path='/projects/ps-glasslab-bioinformatics/software/bowtie2/indexes/'
star_path='/projects/ps-glasslab-bioinformatics/software/STAR/'
bowtie_path='/projects/glass-group/bioinformatics/bowtie2'
homer_path='/projects/glass-group/bioinformatics/homer/bin/'
flip=true

# check number of arguments
if [ $# -lt 4 ] 
then
    echo "Usage: "
    echo "map_on_oasis.sh <experiment type (atac|chip|rna)> <genome> \
<email> <input file directory> [optional arguments]"
    echo "Options:
-l    map files already on tscc or are already copied over
-t    generate qsub scripts but do not execute them
-m    only map files - do not create tag directories
-e    do not send email notifications
-s    copy sam files to Glassome
-p    input data is paired end
-f    do not use the -flip option when making RNA tag directories"
    exit 1
fi

### parse the input ###

OPTIND=5
while getopts "ltmepsf" option ; do # set $o to the next passed option
    case "$option" in  
    l)  
       map_local_files=true 
    ;;  
    t)  
        testing=true
    ;;  
    m)  
        map_only=true
    ;;  
    e)  
        no_emails=true
    ;;  
    p)  
        paired=true
    ;;  
    s)  
        copy_sam=true
    ;;  
    f)  
        flip=false
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
if $paired
then
    echo "Paired end option specified. This script is designed to work with
Illumina paired end reads only"
fi

if $no_emails
then
    echo "Email notifications have been disabled"
else
    echo "Email notifications will be sent to $email"
fi

if $map_only
then
    echo "You are using the map only option - tag directories won't be created
 but sam files will be created on TSCC. Specify the -s option to copy sam files
to Glassome"
fi

if $copy_sam
then
    echo "Sam files will be copied to Glassome"
fi

if $testing
then
    echo "Testing option enabled - qsub scripts will not be activated"
fi

### check arguments ###

#if [ $fileSource == "glassome" ]
if ! $map_local_files
then
    inputDirectory=${inputDirectory#*gpfs/data01/glasslab/}
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
if [ ! $experimentType == "rna" ] && [ ! $experimentType == "chip" ] && [ ! $experimentType == "atac" ]
then
    echo "Error! valid choices for experiment type are 'atac', 'chip' or rna'"
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
#glassomeOutputDirectory="/projects/ps-glasslab-data/scratch/$USER/${inputDirectory##*/}"
glassomeOutputDirectory="$scratch_path/$USER/${inputDirectory##*/}"
if [ ! -d $glassomeOutputDirectory ]
then
    mkdir -p $glassomeOutputDirectory
else
    read -p "This script will copy output files to $glassomeOutputDirectory,
which already exists! Would you like to delete the existing tag directories,
log files, and sam files located at $glassomeOutputDirectory?
Enter y for yes and n for no [yn]" -n 1 -r 
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        rm -rf $glassomeOutputDirectory/log_files/
        rm -rf $glassomeOutputDirectory/sam_files/
        rm -rf $glassomeOutputDirectory/tag_directories
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

# create separate directory for each sample
if [ $(ls $outputDirectory/*fastq.gz| wc -l) -gt 0 ]
then
    for f in $outputDirectory/*fastq.gz;
    do
        dirname=${f%_S[0-9][0-9]*_L0*.fastq.gz}
        dirname=${dirname%_S[0-9][0-9]*.fastq.gz}
        dirname=${dirname%_S[0-9]*.fastq.gz}
        dirname=${dirname%.fastq.gz} # for older data without lane number
        if [ ! -d $dirname ]
        then
            mkdir $dirname
        fi
        mv $f $dirname
    done
fi


if [ $(ls $outputDirectory/*fastq| wc -l) -gt 0 ]
then
    for f in $outputDirectory/*fastq;
    do
        dirname=${f%_S[0-9][0-9]*_L0*.fastq}
        dirname=${dirname%_S[0-9][0-9]*.fastq}
        dirname=${dirname%_S[0-9]*.fastq}
        dirname=${dirname%.fastq} # for older data without lane number
        if [ ! -d $dirname ]
        then
            mkdir $dirname
        fi
        mv $f $dirname
    done
fi


if [ $(ls $outputDirectory/*sra| wc -l) -gt 0 ]
then
    for f in $outputDirectory/*sra;
    do
        dirname=${f%_S[0-9][0-9]*_L0*.sra}
        dirname=${dirname%_S[0-9][0-9]*.sra}
        dirname=${dirname%_S[0-9]*.sra}
        dirname=${dirname%.sra} # for older data without lane number
        if [ ! -d $dirname ]
        then
            mkdir $dirname
        fi
        mv $f $dirname
    done
fi

### decompress files

echo "Decompressing raw data (fastq.gz files)"

# find fastq.gz files
compressedDirs=()
echo $outputDirectory
compressedPaths_1=( $(find $outputDirectory -path "*fastq.gz" -type f) )
compressedPaths_2=( $(find $outputDirectory -path "*sra" -type f) )
compressedPaths_3=( $(find $outputDirectory -path "*fastq" -type f) )
compressedPaths=( ${compressedPaths_1[@]} ${compressedPaths_2[@]} ${compressedPaths_3[@]})
for f in ${compressedPaths[*]}
do
    # remove file name to get sample directory
    compressedDir=${f%/*gz}
    compressedDir=${compressedDir%/*sra}
    compressedDir=${compressedDir%/*fastq}
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
    if $paired
    then
        # for paired end sequence
        if ls $sample_dir/*.fastq &> /dev/null; then
            if [ `ls -l $sample_dir/*.fastq | wc -l` -ne 2 ]; then
                # for first set of reads
                cat $sample_dir/*R1*.fastq > $sample_dir/$bname.fastq_joined
                rm $sample_dir/*R1*.fastq
                mv $sample_dir/$bname.fastq_joined $sample_dir/${bname}_R1.fastq
                # for second set of reads
                cat $sample_dir/*R2*.fastq > $sample_dir/$bname.fastq_joined
                rm $sample_dir/*R2*.fastq
                mv $sample_dir/$bname.fastq_joined $sample_dir/${bname}_R2.fastq
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
                if [ `ls -l $sample_dir/*.fastq | wc -l` -ne 2 ]; then
                    # for first set of reads
                    cat $sample_dir/*R1*.fastq > $sample_dir/$bname.fastq_joined
                    rm $sample_dir/*R1*.fastq
                    mv $sample_dir/$bname.fastq_joined $sample_dir/${bname}_R1.fastq
                    # for second set of reads
                    cat $sample_dir/*R2*.fastq > $sample_dir/$bname.fastq_joined
                    rm $sample_dir/*R2*.fastq
                    mv $sample_dir/$bname.fastq_joined $sample_dir/${bname}_R2.fastq
                else
                    # Rename singular dumped sra file
                    mv $sample_dir/*R1*.fastq $sample_dir/${bname}_R1.fastq
                    mv $sample_dir/*R2*.fastq $sample_dir/${bname}_R2.fastq
                fi  
            fi  
            # Make single file, unzipping simultaneously if they are zipped
            if ls $sample_dir/*.gz &> /dev/null; then 
                # for first set of reads
                zcat $sample_dir/*R1*.gz > $sample_dir/${bname}_R1.fastq
                # for second set of reads
                zcat $sample_dir/*R2*.gz > $sample_dir/${bname}_R2.fastq
            fi  
        fi  
    else
        # for single end sequencing
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
    if [ $(ls $outputDirectory/qsub_scripts/* |wc -l) -ne 0 ]
    then
        rm $outputDirectory/qsub_scripts/*
    fi
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
#if [ ! -d $outputDirectory/pbc ]
#then
#    mkdir $outputDirectory/pbc
#fi

### generate qsub scripts ###

# generate a UUID for this set of jobs
uuid=$[ 1 + $[ RANDOM % 10000 ]] # generate a random number between 0 and 10000
uuid=${USER}_${uuid}

# find directory where script is located
codebase=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) 

# generate script to map fastq to genome

# find all fastq files
echo "generating scripts"
for sampleDir in ${sampleDirs[*]}
    do
    if $paired
    then
        fastqFiles="$(readlink -fm $sampleDir/*R1*fastq) $(readlink -fm $sampleDir/*R2*fastq)"
        fastqFile=$(readlink -fm $sampleDir/*R1*fastq)
        currentDirectory=${fastqFile%/*R1*fastq}
        sampleName=${fastqFile%_R1.fastq}
    else
        fastqFile=$(readlink -fm $sampleDir/*fastq)
        fastqFiles=$fastqFile
        currentDirectory=${fastqFile%/*fastq}
        sampleName=${fastqFile%.fastq}
    fi
    sampleName=${sampleName##/*/} # remove preceding file path
    sampleName=${sampleName#Sample_} # remove "Sample_" from file names"

    samName=""
    logName=""

    # map file
    ### ChIP-seq ###
    if [ $experimentType == "chip" ]
    then
        samName="${sampleName}.${genome}.bowtie2.sam" # change extension to sam
        logName="${sampleName}.${genome}.bowtie2.log" # remove path preceding file name
        if [ -f $outputDirectory/log_files/$logName ]
        then
            rm $outputDirectory/log_files/$logName
        fi

        # execute bowtie
        if $paired
        then
            fastqFiles="-1 ${fastqFiles/ / -2 }"
        fi
        command="$bowtie_path/bowtie2 \
-p 8 \
-x $bowtie_index_path/$genome \
$fastqFiles \
> $outputDirectory/sam_files/$samName \
2> $outputDirectory/log_files/$logName \n"

        # create tag directory
        if ! $map_only
        then
            command+="$homer_path/makeTagDirectory \
$outputDirectory/tag_directories/$sampleName \
-genome $genome \
-checkGC $outputDirectory/sam_files/$samName "
            if $paired
            then
                command+="-sspe "
            fi
            command+="-format sam \n"
        fi

    ### RNA-seq ###
    elif [ $experimentType == "rna" ]
    then
        samName="${sampleName}.${genome}.star.sam" # change extension to sam
        logName="${sampleName}.${genome}.star.log" # remove path preceding file name
        if [ -f $outputDirectory/log_files/$logName ]
        then
            rm $outputDirectory/log_files/$logName
        fi
        # execute star
        command="$star_path/STAR \
--genomeDir $star_path/genomes/$genome \
--readFilesIn $fastqFiles \
--outFileNamePrefix $currentDirectory/ \
--runThreadN 8\n"
        # rename aligned file
        command+="mv $currentDirectory/Aligned.out.sam \
$outputDirectory/sam_files/$samName\n"
        # rename log file
        command+="mv $currentDirectory/Log.final.out \
$outputDirectory/log_files/$logName\n"
        # create tag directory
        if ! $map_only
        then
            command+="$homer_path/makeTagDirectory \
$outputDirectory/tag_directories/${sampleName} \
-genome $genome \
-checkGC $outputDirectory/sam_files/$samName "
            if $paired
            then
                command+="-sspe "
            fi
            if $flip         
            then
                command+="-format sam -flip\n"
            else
                command+="-format sam\n"
            fi
        fi
    ### ATAC-seq ###
    elif [ $experimentType == "atac" ]
    then
        samName="${sampleName}.${genome}.bowtie2.sam" # change extension to sam
        logName="${sampleName}.${genome}.bowtie2.log" # remove path preceding file name
        if [ -f $outputDirectory/log_files/$logName ]
        then
            rm $outputDirectory/log_files/$logName
        fi
        # execute bowtie
        if $paired
        then
            fastqFiles="-1 ${fastqFiles/ / -2 }"
        fi
        command="$bowtie_path/bowtie2 \
-p 8 \
-x $bowtie_index_path/$genome \
$fastqFiles \
> $outputDirectory/sam_files/$samName \
2> $outputDirectory/log_files/$logName \n"
        # create tag directory
        if ! $map_only
        then
        command+="$homer_path/makeTagDirectory \
$outputDirectory/tag_directories/${sampleName}_with_M \
-genome $genome \
-checkGC $outputDirectory/sam_files/$samName "
        if $paired
        then
            command+="-sspe "
        fi
        command+="-format sam\n"
        # remove contaminating tags from chromosome M
        command+="rm $outputDirectory/tag_directories/${sampleName}_with_M/chrM.tags.tsv\n"
        # remake tag directory
        command+="$homer_path/makeTagDirectory \
$outputDirectory/tag_directories/${sampleName} -d \
$outputDirectory/tag_directories/${sampleName}_with_M\n"
        # copy original tag info file
        command+="mv \
$outputDirectory/tag_directories/${sampleName}_with_M/tagInfo.txt \
$outputDirectory/tag_directories/${sampleName}/tagInfo_with_M.txt\n"
        # remove original tag directory
        command+="rm -rf $outputDirectory/tag_directories/${sampleName}_with_M\n"
        fi

    else
        echo "Error! valid choices for experiment type are atac, chip or rna"
        exit 1
    fi

    # calculate PBC coefficient

#    if [ "$experimentType" == "atac" ] || [ "$experimentType" == "chip" ]
#    then
#        uniqueFile=$outputDirectory/pbc/${sampleName}.unique.bam
#        sortedFile=$outputDirectory/pbc/${sampleName}.sorted
#        pileupFile=$outputDirectory/pbc/${sampleName}.pileup
#
#        command+="samtools view -Sbq 1 $outputDirectory/sam_files/$samName > \
#    $uniqueFile\n"
#        command+="samtools sort $uniqueFile $sortedFile\n"
#        command+="samtools mpileup ${sortedFile}.bam > $pileupFile\n"
#        command+="PBC=\$(awk 'BEGIN {N1=0;ND=0} {if(\$4==1){N1+=1} ND+=1} END{print N1/ND}' ${pileupFile})\n"
#        command+="echo -e \"PBC    \$PBC\" >>$outputDirectory/log_files/$logName\n" #"
#    fi

    
    # copy files to Glassome scratch directory
    # copy log file
    command+="cp $outputDirectory/log_files/$logName \
$glassomeOutputDirectory/log_files/\n"
    if ! $map_only
    then
        # copy log file to tag directory
        command+="cp $outputDirectory/log_files/$logName \
$outputDirectory/tag_directories/$sampleName/\n"
        command+="cp -r $outputDirectory/tag_directories/$sampleName \
$glassomeOutputDirectory/tag_directories\n"
    fi

    if $copy_sam
    then
        command+="cp $outputDirectory/sam_files/$samName \
${glassomeOutputDirectory}\n"
    fi
    command+="chmod -R g+w $glassomeOutputDirectory/tag_directories/$sampleName\n"

    # create qsub script
    if $no_emails
    then
    echo -e "#!/bin/bash
#PBS -q hotel
#PBS -N ${sampleName}
#PBS -l nodes=1:ppn=8
#PBS -l walltime=2:00:00
#PBS -o $outputDirectory/qsub_scripts/${sampleName}_torque_output.txt
#PBS -e $outputDirectory/qsub_scripts/${sampleName}_torque_error.txt
#PBS -m n 
#PBS -V
#PBS -A glass-group
$command" > $outputDirectory/qsub_scripts/${sampleName}.torque.sh
    else
    echo -e "#!/bin/bash
#PBS -q hotel
#PBS -N ${sampleName}
#PBS -l nodes=1:ppn=8
#PBS -l walltime=2:00:00
#PBS -o $outputDirectory/qsub_scripts/${sampleName}_torque_output.txt
#PBS -e $outputDirectory/qsub_scripts/${sampleName}_torque_error.txt
#PBS -M $email
#PBS -V
#PBS -m abe
#PBS -A glass-group
$command" > $outputDirectory/qsub_scripts/${sampleName}.torque.sh
    fi
    

    # submit script
    if ! $testing
        then
        echo "Submitting job for $sampleName"
        qsub $outputDirectory/qsub_scripts/${sampleName}.torque.sh
        chmod a+r $outputDirectory/qsub_scripts/*
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

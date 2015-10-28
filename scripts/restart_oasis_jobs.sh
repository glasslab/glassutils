#! /bin/bash
################################################################################

map_local_files=false
testing=false
map_only=false
no_emails=false
paired=false
glassome_path='/projects/ps-glasslab-data/'
mappingScripts_path='/projects/ps-glasslab-bioinformatics/glassutils/mapping_scripts/'
bowtie_index_path='/projects/ps-glasslab-bioinformatics/software/bowtie2/indexes/'
star_path='/projects/ps-glasslab-bioinformatics/software/STAR/'
bowtie_path='/projects/glass-group/bioinformatics/bowtie2'
homer_path='/projects/glass-group/bioinformatics/homer/bin/'
#homer_path='/projects/ps-glasslab-bioinformatics/homer/bin/'

# check number of arguments
if [ $# -lt 1 ] 
then
    echo "Usage: "
    echo "restart_oasis_jobs.sh <path to project directory>"
    exit 1
fi

### parse the input ###


directory=$1

if [ ! -d $directory ]
then
    # directory does not exist - must be a directory on Glassome
    inputDirectory=$(readlink -fm ${glassome_path}/${directory/data//})
    directory="/oasis/tscc/scratch/$USER/${inputDirectory##*/}"
fi

for tagDir in $directory/tag_directories/*
do
    if [ $(ls $tagDir|wc -l) -lt 20 ]
    then
        if [ ! -f $tagDir/tagInfo.txt ]
        then
            qsub $directory/qsub_scripts/${tagDir##*/}.torque.sh
        fi
    fi
done

#!/bin/bash

CMD=$1
GET_DIR=$2
TO_DIR=$3
GENOME=$4
BIOWHAT_USER=$5
EMAIL=$6

EXEC_DIR=/projects/glass-group/glassutils/glassutils/mapping_scripts
BOWTIE_INDEXES=/projects/ps-glasslab/bioinformatics/software/bowtie2/indexes
GTF_FILES=/projects/ps-glasslab/bioinformatics/software/tophat2/gtf

if [ "$GENOME" == "" ]; then
    GENOME="mm9"
fi
if [ "$BIOWHAT_USER" == "" ]; then
    BIOWHAT_USER="$USER"
fi
if [ "$EMAIL" == "" ]; then
    EMAIL="$USER@ucsd.edu"
fi

if [ "$CMD" == "" ] || [ "$CMD" == "help" ]; then
    curr_file=`basename ${0}`
    
    echo "
Please enter a command: send, clean, tophat, bowtie.

Usage: ${curr_file} <1: tophat|bowtie|send|clean > <2: path to data on glassome > <3: where to put data on TSCC > <4: genome > <5: username on glassome > <6: email address >

For help, please see https://sites.google.com/a/glasso.me/glasslab/computational-guides/using-map_on_tscc-sh"
    exit
fi


# Get File locations, using readlink to get rid of double-slashes,
# since those break fastq-dump
DATA_DIR=$(readlink -m $TO_DIR/`basename $GET_DIR`)

if [ "$CMD" == "send" ]; then
    echo "Copying ${DATA_DIR} to ${BIOWHAT_USER}@glassome.ucsd.edu:${GET_DIR}..."
    mkdir $DATA_DIR/processed
    mv $DATA_DIR/*/*${GENOME}* $DATA_DIR/processed
    cp -r $DATA_DIR/processed /projects/ps-glasslab/$GET_DIR 
    exit
else
    if [ "$CMD" == "clean" ]; then
        if [ "$DATA_DIR" == "" ]; then
            echo "Cannot delete specified directory."
            exit
        fi
        echo -n "Are you sure you want to remove the entire directory ${DATA_DIR}? [Y/n] "
        read -e CONFIRM
        if [ "$CONFIRM" == "Y" ] || [ "$CONFIRM" == "y" ] || [ "$CONFIRM" == "yes" ]; then
            echo "Deleting ${DATA_DIR}."
            rm -r $DATA_DIR
            exit
        else
            echo "Not deleting ${DATA_DIR}. Exiting."
            exit
        fi
    else
        # Set up operation
        if [ "$CMD" == "tophat" ]; then
            OP="perl ${EXEC_DIR}/map-bowtie2.pl -index ${BOWTIE_INDEXES}/${GENOME} -cpu 1 -p 8 --library-type fr-secondstrand -G ${GTF_FILES}/${GENOME}.refseq.gtf -tophat2 "
            WTIME="40:00:00"
            NODES="nodes=1:ppn=8"
        else
            if [ "$CMD" == "bowtie" ]; then
                OP="perl ${EXEC_DIR}/map-bowtie2.pl -index ${BOWTIE_INDEXES}/${GENOME} -cpu 1 -p 8 "
                WTIME="3:00:00"
                NODES="nodes=1:ppn=8"
            else
		if [ "$CMD" == "STAR" ]; then
			OP="perl ${EXEC_DIR}/map-star.pl ${GENOME} "
			WTIME="3:00:00"
			NODES="nodes=1:ppn=8"
		else
			echo "Did not recognize command '${CMD}'. Exiting."
			exit
		fi
            fi
        fi        
        # Move files over
        if [ "${GET_DIR:0:5}" != "local" ]; then
            cp -r /projects/ps-glasslab/$GET_DIR $TO_DIR 
        fi


        # Decompress files, combine.
        for SUB_DIR in $DATA_DIR/*
          do
            bname=`basename $SUB_DIR`

            # If there is exactly one fastq file, use that; otherwise...
            if ls $SUB_DIR/*.fastq &> /dev/null; then
                if [ `ls -l $SUB_DIR/*.fastq | wc -l` -ne 1 ]; then
                    cat $SUB_DIR/*.fastq > $SUB_DIR/$bname.fastq_joined
                    rm $SUB_DIR/*.fastq
                    mv $SUB_DIR/$bname.fastq_joined $SUB_DIR/$bname.fastq
                fi
                # else, only one .fastq; will be used.
            else
                # If there are any .sra files, dump to .fastq
                if ls $SUB_DIR/*.sra &> /dev/null; then 
                    for sra in $SUB_DIR/*.sra
                        do
                            CURR_DIR=`pwd`
                            # CD in so that fastq-dump works correctly
                            cd $SUB_DIR
                            fastq-dump $sra
                            rm $sra
                            cd $CURR_DIR
                        done
                    # Then compile all the .fastq
                    if [ `ls -l $SUB_DIR/*.fastq | wc -l` -ne 1 ]; then
                        cat $SUB_DIR/*.fastq > $SUB_DIR/$bname.fastq_joined
                        rm $SUB_DIR/*.fastq
                        mv $SUB_DIR/$bname.fastq_joined $SUB_DIR/$bname.fastq
                    else
                        # Rename singular dumped sra file
                        mv $SUB_DIR/*.fastq $SUB_DIR/$bname.fastq
                    fi
                fi
                # Make single file, unzipping simultaneously if they are zipped
                if ls $SUB_DIR/*.gz &> /dev/null; then 
                    zcat $SUB_DIR/*.gz > $SUB_DIR/$bname.fastq
                fi
            fi
          done

        # Create PBS file for each
        for fastq in $DATA_DIR/*/*.fastq
          do
            OP_for_file="${OP} ${fastq}"
            bname_fastq=`basename $fastq`
            job_file=$DATA_DIR/${bname_fastq}_job_file.sh
	    echo $job_file

if [ "$CMD" == "STAR" ]; then
	new_dir=`dirname $fastq`
else
	new_dir="/oasis/tscc/scratch/${USER}"
fi

# Note that leading whitespace breaks Torque. 
echo "#!/bin/bash
#PBS -q hotel
#PBS -N ${bname_fastq}
#PBS -l ${NODES}
#PBS -l walltime=${WTIME}
#PBS -o ${fastq}.torque_output.txt
#PBS -e ${fastq}.torque_error.txt
#PBS -V
#PBS -M ${EMAIL}
#PBS -m abe
#PBS -A glass-group

cd $new_dir

${OP_for_file}" > $job_file

            qsub $job_file
        done
        exit
    fi
fi

CMD=$1
GET_DIR=$2
TO_DIR=$3
GENOME=$4
BIOWHAT_USER=$5
EMAIL=$6

EXEC_DIR=/projects/glass-lab/bioinformatics
BOWTIE_INDEXES=/projects/glass-lab/bioinformatics/bowtie2/indexes
GTF_FILES=/projects/glass-lab/bioinformatics/tophat2/gtf

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
	echo "Please enter a command: send, clean, tophat, bowtie."
	curr_file = `basename ${0}`
	echo "Usage: ${curr_file} <1: send|clean|tophat|bowtie > <2: path to data on biowhat > <3: where to put data on Triton > <4: genome > <5: username on biowhat > <6: email address >"
	echo "For help, please contact karmelallison@ucsd.edu"
	exit
fi


# Get File locations
DATA_DIR=$TO_DIR/`basename $GET_DIR`

if [ "$CMD" == "send" ]; then
	echo "Copying ${DATA_DIR} to ${BIOWHAT_USER}@biowhat.ucsd.edu:${GET_DIR}..."
	scp -r $DATA_DIR/*/*${GENOME}* $BIOWHAT_USER@biowhat.ucsd.edu:$GET_DIR 
	exit
else
	if [ "$CMD" == "clean" ]; then
		if [ "$DATA_DIR" == "" ]; then
			echo "Cannot delete specified directory."
			exit
		fi
		echo -n "Are you sure you want to remove the entire directory ${DATA_DIR}? [Y/n]"
		read -e CONFIRM
		if [ "$CONFIRM" == "Y" ] || [ "$CONFIRM" == "y" ] || [ "$CONFIRM" == "yes" ]; then
			echo "Deleting ${DATA_DIR}."
			#rm -r $DATA_DIR
			exit
		else
			echo "Not deleting ${DATA_DIR}. Exiting."
			exit
		fi
	else
		# Set up operation
		if [ "$CMD" == "tophat" ]; then
			OP="perl ${EXEC_DIR}/map-bowtie2.pl -index ${BOWTIE_INDEXES}/${GENOME} -cpu 1 -p 1 --library-type fr-secondstrand -G ${GTF_FILES}/${GENOME}.refseq.gtf -tophat2 ${fastq}"
			WTIME="10:00:00"
			NODES="nodes=1:ppn=1"
		else
			if [ "$CMD" == "bowtie" ]; then
				OP="perl ${EXEC_DIR}/map-bowtie2.pl -index ${BOWTIE_INDEXES}/${GENOME} -cpu 1 -p 8 ${fastq}"
				WTIME="3:00:00"
				NODES="nodes=1:ppn=8"
			else
				echo "Did not recognize command '${CMD}'. Exiting."
				exit
			fi
		fi
		
		# Move files over
		scp -r $BIOWHAT_USER@biowhat.ucsd.edu:$GET_DIR $TO_DIR 


		# Decompress files, combine.
		for SUB_DIR in $DATA_DIR/*
		  do
		    bname=`basename $SUB_DIR`
		    zcat $SUB_DIR/*.gz > $SUB_DIR/$bname.fq
		  done

		# Create PBS file for each
		for fastq in $DATA_DIR/*/*.fq
		  do
			bname_fq=`basename $fastq`
			job_file=$DATA_DIR/${bname_fq}_job_file.sh
			echo "#!/bin/bash
			#PBS -q small
			#PBS -N ${bname_fq}
			#PBS -l ${NODES}
			#PBS -l walltime=${WTIME}
			#PBS -o ${fastq}.torque_output.txt
			#PBS -e ${fastq}.torque_error.txt
			#PBS -V
			#PBS -M ${EMAIL}
			#PBS -m abe
			#PBS -A glass-lab

			cd /phase1/${USER}

			${OP}" > $job_file

			qsub $job_file
		done
		exit
	fi
fi
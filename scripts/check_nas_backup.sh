#!/bin/bash

TEST_PATH=/data/timed_test_file_for_backup_monitoring.txt
BACKUP_PATH=/share/backup/hourly.0/biowhat.ucsd.edu

# First, touch target file for next time.
touch $TEST_PATH

# Then confirm that the file on the backup server has been modified in the last day
CMD="find ${BACKUP_PATH}${TEST_PATH} -mtime -1"
CONFIRM=`ssh admin@nas.glasso.me ${CMD}`

if [ "$CONFIRM" != "${BACKUP_PATH}${TEST_PATH}" ]; then
	# Nagios alert
	echo "No"
fi
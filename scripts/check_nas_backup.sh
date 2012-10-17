#!/bin/bash

TEST_PATH=/data/timed_test_file_for_backup_monitoring.txt
BACKUP_PATH=/share/backup/hourly.0/biowhat.ucsd.edu

# First, touch target file for next time.
touch $TEST_PATH

# Then confirm that the file on the backup server has been modified in the last day
CMD="find ${BACKUP_PATH}${TEST_PATH} -mtime -1"
CONFIRM=`ssh admin@nas.glasso.me ${CMD}`

if [ "${CONFIRM}" != "${BACKUP_PATH}${TEST_PATH}" ]; then
	# Nagios alert- critical
	echo "Did not recognize NAS backup test file.
Ran command: $CMD
Found: ${CONFIRM}
Was looking for: ${BACKUP_PATH}${TEST_PATH}"
	exit 2
fi

echo "Found NAS backup test file: ${CONFIRM} == ${BACKUP_PATH}${TEST_PATH}"
exit 0
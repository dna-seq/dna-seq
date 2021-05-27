#!/bin/bash
ABSPATH=$(readlink -f $0)
ABSDIR=$(dirname $ABSPATH)
if ! [ $(id -u) = 0 ]; then
	echo "Must be root!"
	exit 1
fi
if ! [ "$1" == "all" ]; then
	echo "Dedicated VMs only!"
	echo "Will wipe docker and cromwell data to the clean slate!"
	echo "usage: purge.sh all"
	exit 1
fi
cd ${ABSDIR}/../
docker system prune -a
rm -rf ${ABSDIR}/../data/cromwell-executions/
rm -rf ${ABSDIR}/../data/cromwell-workflow-logs/workflow.*
rm -rf ${ABSDIR}/../data/databases/mysql/*
git restore ./data/cromwell-executions/.gitignore
git restore ./data/cromwell-workflow-logs/.gitignore
git restore ./data/cromwell-workflow-logs/.gitignore
cd -

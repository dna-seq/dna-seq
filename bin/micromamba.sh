#!/bin/bash
ABSPATH=$(readlink -f $0)
ABSDIR=$(dirname $ABSPATH)
if [ $(id -u) = 0 ]; then
   echo "Must be user!"
   exit 1
fi

if ! [ -d "${ABSDIR}/../data/micromamba" ]; then
	micromamba shell init -s bash -p "${ABSDIR}/../data/micromamba"
	source ~/.bashrc
	echo $MAMBA_ROOT_PREFIX
	micromamba create -f "${ABSDIR}/../environment.yaml"
fi

echo Complete!

#!/bin/bash
################################ Server Specs
#PBS -j oe
#PBS -m n
#PBS -N classical
#PBS -q serial
#PBS -l nodes=1:ppn=1
################################ Meat & Potatoes
mkdir -p "/scratch/$USER" 			# make sure the scratch directory is available
echo "Running script on $HOSTNAME"	# simple output
/home/ngraymon/stable_mmtk/bin/python "$PBS_O_WORKDIR"/test.py -d "$PBS_O_WORKDIR/workspace/" -f "betaho$NUM"

mv --force "/scratch/$USER/classicalbetaho$NUM"* "$PBS_O_WORKDIR/workspace"
 
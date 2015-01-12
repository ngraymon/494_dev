#!/bin/bash

################################ Server Stuff ################################ 

# this command joins the stderr and stdout
#PBS -j oe

# this command stops all mail from being sent
#PBS -m n

# this command give the job a name
#PBS -N simple_pigs_HO

# this command defines the priority of the job, default is zero
#PBS -p 0

# this command specifies the queue to submit to
#PBS -q long

# this command request # of nodes and # of processor per node
#PBS -l nodes=1:ppn=1

################################ Meat & Potatoes ################################ 

# make sure the scratch directory is available
mkdir -p "/scratch/$USER"

echo "Running script on $HOSTNAME"
echo


if   [[ -n "$BETA" ]]; then
	/home/ngraymon/stable_mmtk_pigs/bin/python "$PBS_O_WORKDIR"/pathbuilder_PIGS.py -d "/scratch/$USER/" -f "beta$FILENAME" -b "$BETA"
	echo "Output file is here: $HOSTNAME:/scratch/$USER/beta$FILENAME.nc"
	mv --force "/scratch/$USER/beta$FILENAME"* "$PBS_O_WORKDIR/workspace"
elif [[ -n "$TAU"  ]]; then
	/home/ngraymon/stable_mmtk_pigs/bin/python "$PBS_O_WORKDIR"/pathbuilder_PIGS.py -d "/scratch/$USER/" -f "tau$FILENAME" -t "$TAU"
	echo "Output file is here: $HOSTNAME:/scratch/$USER/tau$FILENAME.nc"
	mv --force "/scratch/$USER/tau$FILENAME"* "$PBS_O_WORKDIR/workspace"
elif [[ -n "$EST"  ]]; then
	/home/ngraymon/stable_mmtk_pigs/bin/python "$PBS_O_WORKDIR"/estimator.py -d "/scratch/$USER/" -f "beta$FILENAME" -b "$BETA"
	echo "Output file is here: $HOSTNAME:/scratch/$USER/beta$FILENAME.plt"
	mv --force "/scratch/$USER/beta$FILENAME"* "$PBS_O_WORKDIR/workspace"
else
	/home/ngraymon/stable_mmtk_pigs/bin/python "$PBS_O_WORKDIR"/pathbuilder_PIGS.py -d "/scratch/$USER/" -f "$FILENAME"
	echo "Output file is here: $HOSTNAME:/scratch/$USER/$FILENAME.nc"
	mv --force "/scratch/$USER/FILENAME"* "$PBS_O_WORKDIR/workspace"
fi


# and done
cd "$PBS_O_WORKDIR"
echo Job Script Finished Successfully 
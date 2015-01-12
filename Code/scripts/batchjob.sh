#!/bin/bash

################################################################################################################
if [[ "$#" -eq 0 || "$1" == "-h" ]]; then
	echo "Possible args :"
	echo "                1) Provide a single integer argument that is the number of jobs you wish to submit"
	echo "                2) Provide 'tau' and two integer arguments that are the limits on a sequence of trial 'tau' values"
	echo "                2) Provide 'tau' and three integer arguments start, finish, index, for a sequence of trial 'tau' values"
	echo "                3) Provide 'beta' and two integer arguments that are the limits on a sequence of trial 'beta' values"
	echo "                3) Provide 'beta' and three integer arguments start, finish, index, for a sequence of trial 'beta' values"
	exit
fi

if [[ "$#" -eq 1 ]]; then
	for (( num = 1; num <=${1}; num++ )) do
		qsub -N "ho$num" -v "FILENAME=ho$num" submit_script.sh
	done
	exit
fi

################################################################################################################

if   [[ "$#" -eq 3  &&  "$1" == "tau" ]]; then
	for (( num = ${2}; num <=${3}; num++ ))  do
		qsub -N "tauho$num" -v "FILENAME=ho$num","TAU=$num" submit_script.sh
	done
	printf '%s\n%s\n' $2 $3 $1 > plotdata.batchjoboutput
	exit
elif [[ "$#" -eq 3  &&  "$1" == "beta" ]]; then
	for (( num = ${2}; num <=${3}; num++ )) do
		qsub -N "betaho$num" -v "FILENAME=ho$num","BETA=$num" submit_script.sh
	done
	printf '%s\n%s\n' $2 $3 $1 > plotdata.batchjoboutput
	exit
fi

################################################################################################################
if [[ "$3" < 1 ]]; then
	if   [[ "$#" -eq 4  &&  "$1" == "tau" ]]; then
		for (( i = $(bc<<<"$2/$4"); i <=$(bc<<<"$3/$4"); i ++ )) do
			NUM=$(bc<<<"$4 * $i")
			qsub -N "tauho$NUM" -v "FILENAME=ho$NUM","TAU=$NUM" submit_script.sh
		done
		printf '%s\n%s\n%s\n' $2 $3 $4 $1 > plotdata.batchjoboutput
		exit
	elif [[ "$#" -eq 4  &&  "$1" == "beta" ]]; then
		for (( i = $(bc<<<"$2/$4"); i <=$(bc<<<"$3/$4"); i ++ )) do
			NUM=$(bc<<<"$4 * $i")
			qsub -N "betaho$NUM" -v "FILENAME=ho$NUM","BETA=$NUM" submit_script.sh
		done
		printf '%s\n%s\n%s\n' $2 $3 $4 $1 > plotdata.batchjoboutput
		exit
	fi

##############

elif [[ ! "$3" < 1 ]]; then
	if   [[ "$#" -eq 4  &&  "$1" == "tau" ]]; then
		for (( num = $2; num <=$3; num += $4 )) do
		qsub -N "tauho$num" -v "FILENAME=ho$num","TAU=$num" submit_script.sh
		done
		printf '%s\n%s\n%s\n' $2 $3 $4 $1 > plotdata.batchjoboutput
		exit
	elif [[ "$#" -eq 4  &&  "$1" == "beta" ]]; then
		for (( num = $2; num <=$3; num += $4 )) do
		qsub -N "betaho$num" -v "FILENAME=ho$num","BETA=$num" submit_script.sh
		done
		printf '%s\n%s\n%s\n' $2 $3 $4 $1 > plotdata.batchjoboutput
		exit
	fi
fi

################################################################################################################

if [[ "$1" != "tau" && "$1" != "beta" ]]; then
	echo "You forgot to provide 'tau' or 'beta' as the first argument"
else
	echo "Incorrect args"
fi

exit
#!/bin/bash
# Submit from a clean environment!
# The scripts set up the environment!
#R_ANA_NAME='SuperflowQFlip'
#R_ANA_NAME='SuperflowFakes'
R_ANA_NAME='runTruthAna'
R_START_DIR='/gdata/atlas/dantrim/SusyAna/Super/'
R_OUTPUT_DIR='/gdata/atlas/dantrim/SusyAna/histoAna/TAnaOutput/SMCwslep/Truth_Apr27/Raw/'
R_LOG_DIR='/gdata/atlas/dantrim/SusyAna/histoAna/TAnaOutput/SMCwslep/Truth_Apr27/logs/'

R_WORK_BASE='/scratch/dantrim/'

R_RUN='/gdata/atlas/dantrim/SusyAna/Super/TAnalysis/run/'
R_GEN='/gdata/atlas/dantrim/SusyAna/Super/TAnalysis/run/filelists/'
R_LIST_DIR=${R_GEN}
R_LIST_POST='.txt'
#R_LIST_POST='_n0145.txt'

R_SAMPLES_LIST=("file_list_SMCwSlep_n0155")


for file_ in ${R_SAMPLES_LIST[@]}; do
	FILE_NAME=${R_LIST_DIR}${R_LIST_PREFIX}${file_}${R_LIST_POST}
	FILE_LINES=`cat $FILE_NAME`
	
	for line in $FILE_LINES ; do
        	export S_ANA_NAME=${R_ANA_NAME}
#		export S_MODE='t' # run 2D param
                export S_MODE='c'
		export S_STARTDIR=${R_START_DIR}
		
		sleep 0.13 # to regenerate entropy
		RANNUM=$RANDOM$RANDOM$RANDOM
		
		export S_WORKDIR=${R_WORK_BASE}${RANNUM}
		export S_IN_DIRECTORY=${line%?}'/'
	#	export S_IN_DIRECTORY=${line%?}'t'
		export S_SYSTEMATIC='NONE'
		export S_OUTPUT_DIR=${R_OUTPUT_DIR}
		
		RED_line=${line%?}
		lFileName=$(basename $RED_line)
		strip_one=${lFileName#user.*.}
		strip_two=${strip_one%.SusyNt*}
		
		echo $strip_two
		
		sbatch -J 'truthRazor '${strip_two} -o ${R_LOG_DIR}${lFileName}_slurm-%j.log /gdata/atlas/dantrim/SusyAna/Super/TruthAna/TruthSelector.sh 
	done
done

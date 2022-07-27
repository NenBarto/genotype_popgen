#!/usr/local/bin/bash

#create and start an environment

#conda create --name enviroDNA
#conda activate enviroDNA
##set up channels
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge

#conda install vsearch 
conda activate enviroDNA
COL_RED=""
COL_RESET=""
COL_WHITE=""
COL_YELLOW=""

projectDir="/share/ScratchGeneral/nenbar/projects/Anthony"
inDir=$projectDir/"raw_data/210813_miseq_AR/fastq/fastq"
#first make all the scripts in bash
SCRIPTDIR=$projectDir"/scripts"
source "$SCRIPTDIR"/_ALL_FUNCTIONS.sh --source-only

run_derep="$SCRIPTDIR/derep_fa.sh"
run_vsearch="$SCRIPTDIR/vsearch_pipeline_v5modified.sh"

WORKDIR=$inDir

trimfunc="merged"
amp="mhc"
db="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/contig2.fasta"
minuniq=10
#echo -e $COL_WHITE "WORKDIR=$WORKDIR" $COL_WHITE
#echo -e $COL_YELLOW "Continue? y|n" $COL_RESET
#read check_workdir
#[[ "$check_workdir" == "y" ]] || exit
#
#echo -n "Amplicon? "
#read amp
#echo -n "Trimming merged|single|1step? "
#read trimfunc
#echo -n "Minimum OTU size? "
#read minuniq
#echo -n "Iteration? "
#read i

mkdir -p "$WORKDIR"/unzip_files
cp $WORKDIR/*.gz "$WORKDIR"/unzip_files
gunzip "$WORKDIR"/unzip_files/*

#unzip_files "$WORKDIR"/Fastq fastq gz "$WORKDIR"/unzip_files 2>&1 | tee "$WORKDIR"/unzip.log
print_valid_samples "$WORKDIR"/unzip_files fastq "$WORKDIR"/valid_samples_1.txt

filter_fq "$WORKDIR"/unzip_files "$WORKDIR"/valid_samples_1.txt "$WORKDIR"/filtered_fq 2>&1 | tee "$WORKDIR"/filter_fq.log
print_valid_samples "$WORKDIR"/filtered_fq fastq "$WORKDIR"/valid_samples_2.txt


merge_fq "$WORKDIR"/unzip_files "$WORKDIR"/valid_samples_1.txt "$WORKDIR"/merged_fq 2>&1 | tee "$WORKDIR"/merge_fq.log	
trim_primer "$WORKDIR"/merged_fq "$WORKDIR"/trim_primer "$amp" "$trimfunc" fasta 2>&1 | tee "$WORKDIR"/trim_primer_"$amp"_run"$i".log

#fix here with the new script
"$run_derep" "$WORKDIR" "$amp" "$WORKDIR"/trim_primer_"$amp" "$minuniq" 2>&1| tee "$WORKDIR"/derep_fa"$amp"_run"$i".log

#fix here with a new script 
"$run_vsearch" \
"$WORKDIR" \
"$WORKDIR"/derep_fa"$amp"/seqs.fna \
"$db" "$minuniq" 2>&1| tee "$WORKDIR"/vsearch_"$amp"_run"$i".log


#issues:
#1. the clustering is done on wrongly named files?
#head ../raw_data/210813_miseq_AR/fastq/fastq/size10-sintax95/clustering.uc

#2. the output file is weird
#head ../raw_data/210813_miseq_AR/fastq/fastq/size10-sintax95/all.otutab.txt


#3. the derep files start as they should
#head ../raw_data/210813_miseq_AR/fastq/fastq/derep_famhc/F0933_S314_L001_mergedR1R2_trimmed.derep









#then turn individual tools into containers

#then build them into nextflow

#then transfer to aws



#fastq_merge.sh
#step 1. vsearch mergepairs

#step 2. vserach fastx_revcomp

#step 3. vsearch fastq_join

#step 4. fastqCombinePairedEnd.py


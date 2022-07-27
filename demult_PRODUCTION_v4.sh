#!/usr/local/bin/bash

#Workflow for demultiplexed files from iSeq

#Variables for iteration only. Individual parameters for VSEARCH still need to be set from within script.

SCRIPTDIR="/Users/svs/CESAR-EnviroDNA/scripts"
source "$SCRIPTDIR"/_ALL_FUNCTIONS.sh --source-only

run_derep="$SCRIPTDIR/ngs/derep_fa.sh"
run_vsearch="$SCRIPTDIR/ngs/vsearch_pipeline_v5TEST.sh"


#WORKDIR="/Users/svs/Downloads/Ant-MPP/COI_20210603_074757"
WORKDIR="/Users/svs/Downloads/SmallBatch"


echo -e $COL_WHITE "WORKDIR=$WORKDIR" $COL_WHITE
echo -e $COL_YELLOW "Continue? y|n" $COL_RESET
read check_workdir
[[ "$check_workdir" == "y" ]] || exit


print_valid_amplicons

echo -e $COL_YELLOW "Please provide input for the following questions:" $COL_RESET

echo -n "Amplicon? "
read amp
echo -n "Trimming merged|single|1step? "
read trimfunc
echo -n "Minimum OTU size? "
read minuniq
echo -n "Iteration? "
read i

case "$amp" in

	amphibian)
		db="/Users/svs/Downloads/Wetlands/Wetland-AMPH/REFDB/refseq_AMPH_full_multsp_v3.fasta"
		;;
	
	bryobia-28s)
		db="/Users/svs/Downloads/Bryobia/refdb/Bryobia_28S_ref.fasta"
		;;
	
	bryobia-coi-f1*)
		db="/Users/svs/Downloads/200602_M04416_0083_000000000-CGPNL/BRY/COI_DummyRef.fasta"
		;;
		
	coi-zbj)
		db="/Users/svs/Downloads/Ant-MPP/REFDBs/COIZ_refseq_v3.fasta"
		#db="/Users/svs/Downloads/SpectrumCOIZ/refseqs/refseq_USETHIS_v3.uniq.fasta"
		;;
		
	decar1)
		db="/Users/svs/Downloads/TESTING/DB/UPDATED_SINTAX_FILES/make_sintax_taxassigned/DecaR2_distinct.sintaxXX.fasta.sintaxFULL.fasta"
		;;
	
	decar2)
		db="/Users/svs/Downloads/TESTING/DB/UPDATED_SINTAX_FILES/make_sintax_taxassigned/DecaR2_distinct.sintaxXX.fasta.sintaxFULL.fasta"
		;;
	
	human-hvr2)
		db=""
		;;

	invert-16s)
		#db="/Users/svs/Downloads/TESTING/InvertDB/make_invert3/16S_RANKEDtax2seq_v1.fasta"
		db="/Users/svs/Downloads/MonbulkHorse16S_2021/REFDB/refdb_v3c_unassig.fasta"
		;;
		
	its2)
		db="/Users/svs/Downloads/Ant trnL ITS plant pilot/ITS_refsequences1.trimmed.subsetREF.fasta"
		#db="/Users/svs/Downloads/Ant trnL ITS plant pilot/ITS2.uniq.fasta"
		;;
		
	kdr-2f2r)
		db="/Users/svs/Downloads/RLEM-kdr/Dummy_RefDB.fasta"	
		;;
			
	mhc)
		db="/Users/svs/Downloads/Ant MHC STF/refcontigs.fasta"
		;;
		
	mideca)
		db="/Users/svs/Downloads/TESTING/MiDeca/MiDeca_v1.trimmed.RANKED.fasta"
		;;
	
	mifish-*)
		db="/Users/svs/Downloads/TESTING/DB/UPDATED_SINTAX_FILES/make_sintax_taxassigned/MiFish_distinct.sintaxXX.fasta.sintaxFULL.fasta"
		;;

	teleo)
		db="/Users/svs/Downloads/TESTING/DB/UPDATED_SINTAX_FILES/make_sintax_taxassigned/Teleo_distinct.sintaxXX.fasta.sintaxFULL_v3.fasta"
		;;
	
	vert)
		db="/Users/svs/Downloads/TESTING/DB/UPDATED_SINTAX_FILES/make_sintax_taxassigned/Vert_distinct.sintaxXX.fasta.sintaxFULL.fasta"
		;;
	
	trnl)
		db="/Users/svs/Downloads/Ant-MPP/REFDBs/TrnL_refseq_v1.fasta"
		;;
	
	*)
		echo -e $COL_WHITE "Exiting due to invalid/ ambiguous input." $COL_RESET
		exit
		;;
	
esac



#concat_runs "$WORKDIR"/Fastq-run1 "$WORKDIR"/Fastq-run2 "$WORKDIR"/Fastq 2>&1 | tee "$WORKDIR"/concat.log

unzip_files "$WORKDIR"/Fastq fastq gz "$WORKDIR"/unzip_files 2>&1 | tee "$WORKDIR"/unzip.log
print_valid_samples "$WORKDIR"/unzip_files fastq "$WORKDIR"/valid_samples_1.txt

filter_fq "$WORKDIR"/unzip_files "$WORKDIR"/valid_samples_1.txt "$WORKDIR"/filtered_fq 2>&1 | tee "$WORKDIR"/filter_fq.log
print_valid_samples "$WORKDIR"/filtered_fq fastq "$WORKDIR"/valid_samples_2.txt
#join_fq "$WORKDIR"/filtered_fq "$WORKDIR"/valid_samples_2.txt "$WORKDIR"/joined_fq

case "$trimfunc" in

	single)
		trim_primer "$WORKDIR"/unzip_files "$WORKDIR"/trim_primer "$amp" "$trimfunc" fasta 2>&1 | tee "$WORKDIR"/trim_primer_"$amp"_run"$i".log
		;;
		
	1step)
		merge_fq "$WORKDIR"/unzip_files "$WORKDIR"/valid_samples_1.txt "$WORKDIR"/merged_fq 2>&1 | tee "$WORKDIR"/merge_fq.log
		#print_valid_samples "$WORKDIR"/merged_fq fastq "$WORKDIR"/valid_samples_2.txt
		"/Users/svs/Downloads/Odonata_and_Feral_Horse_VERT_1-Step/demult_1STEP.sh" "$WORKDIR"/merged_fq "$WORKDIR"/demult
		trim_primer "$WORKDIR"/demult "$WORKDIR"/trim_primer "$amp" merged fasta 2>&1 | tee "$WORKDIR"/trim_primer_"$amp"_run"$i".log
		;;
		
	merged)
		merge_fq "$WORKDIR"/unzip_files "$WORKDIR"/valid_samples_1.txt "$WORKDIR"/merged_fq 2>&1 | tee "$WORKDIR"/merge_fq.log	
		trim_primer "$WORKDIR"/merged_fq "$WORKDIR"/trim_primer "$amp" "$trimfunc" fasta 2>&1 | tee "$WORKDIR"/trim_primer_"$amp"_run"$i".log
		;;
		
	*)
		echo -e $COL_WHITE "Exiting due to invalid/ ambiguous input." $COL_RESET
		exit
		;;
esac


#Usage: [WORKDIR] [amp] [trimmed-fasta-dir] (opt:MINUNIQ, def=10)
"$run_derep" "$WORKDIR" "$amp" "$WORKDIR"/trim_primer_"$amp" "$minuniq" 2>&1| tee "$WORKDIR"/derep_fa"$amp"_run"$i".log



#Usage: [WORKDIR] [merged-derep-fasta] [flagged-taxDB-fasta] (opt:SIZE_FILTER, def=10)
#"$run_vsearch" \
#"$WORKDIR" \
#"$WORKDIR"/derep_fa"$amp"/seqs.fna \
#"$db" "$minuniq" 2>&1| tee "$WORKDIR"/vsearch_"$amp"_run"$i".log


get_unassigned_fasta "$WORKDIR"/size10-sintax95/counts_with_tax.txt 50 \
"$WORKDIR"/size10-sintax95/all.otus.fasta \
"$WORKDIR"/Unassigned.fasta

exit


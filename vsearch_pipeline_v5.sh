#!/bin/bash

#v5: added
#--maxaccepts
#--maxrejects
#from previous:
#--clusterout_id
#--userfields

#$() is equivalent to backticks which have now been deprecated

#unofficial bash strict mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'

COL_RED=""


function usage {
    echo "Usage: [WORKDIR] [merged-derep-fasta] [flagged-taxDB-fasta] (opt:size-filter, def=10)"
    exit
}


if [ $# -lt 3 ]
then
    usage
fi


WORKDIR="$1"
projectDir="/share/ScratchGeneral/nenbar/projects/Anthony"
inDir=$projectDir/"raw_data/210813_miseq_AR/fastq/fastq"
#first make all the scripts in 
OUTDIR="$WORKDIR/size10-sintax95"
VSEARCH=$(which vsearch)
THREADS=2


#input
fi_STARTING_FASTA="$2"
fi_TAXDB_FASTA="$3"
SIZE_FILTER="${4:-10}"


#filtering and formatting
fo_STARTING_FASTA="$OUTDIR/starting.fasta"
FASTA_WIDTH=0 #no wrapping

#chimera detection - de novo and ref
fo_ALL_DENOVO_NONCHIMERAS_FASTA="$OUTDIR/all.denovo.nonchimeras.fasta"
fo_ALL_NONCHIMERAS_FASTA="$OUTDIR/all.nonchimeras.fasta"
fo_ALL_NONCHIMERAS_UC="$OUTDIR/all.nonchimeras.uc"

#clustering - options: cluster_size, cluster_unoise, cluster_fast etc.
IDENT=1.0
IDDEF=4
CLUSTER_METHOD="cluster_size"
fo_ALL_OTUS_FASTA="$OUTDIR/all.otus.fasta"
fo_ALL_OTUTAB_TXT="$OUTDIR/all.otutab.txt"
fo_CLUSTERING_UC="$OUTDIR/clustering.uc"

#taxonomic classification - cutoff is min bootstrap support for ranks that will be included
SINTAX_CUTOFF=0.95
fo_ALL_TAXED_OUT="$OUTDIR/all.taxed.out.tsv"
fo_COUNTS_WITH_TAX="$OUTDIR/counts_with_tax.txt"
fo_AGG_COUNTS="$OUTDIR/agg_counts.txt"


if [ ! -d $OUTDIR ]
then
    mkdir $OUTDIR
fi



echo "Input fasta file with sequence-abundance per sample, singletons removed: "$fi_STARTING_FASTA""
echo "Ref: "$fi_TAXDB_FASTA""

echo "No. of sequences in "$fi_STARTING_FASTA": $(grep -c "^>" "$fi_STARTING_FASTA")"


echo
echo "Size (abundance) filtering, minsize=$SIZE_FILTER"

#make single line fasta
sed 's/>/,>/g' "$fi_STARTING_FASTA" | tr '\n' '@' | tr ',' '\n' > "$fo_STARTING_FASTA".tmp

i=$SIZE_FILTER

while [ $i -gt 0 ]
do
	(( i-- ))
	grep -v -w "size=$i" "$fo_STARTING_FASTA".tmp > "$fo_STARTING_FASTA".filtered
	mv "$fo_STARTING_FASTA".filtered "$fo_STARTING_FASTA".tmp
	
done

tr '@' '\n' < "$fo_STARTING_FASTA".tmp | tr -s '\n' '\n' | tr -s '\r' '\n' | tail -n +2 > "$fo_STARTING_FASTA".filtered

echo "Sequences remaining: $(grep -c "^>" "$fo_STARTING_FASTA".filtered)"



cp "$fo_STARTING_FASTA".filtered "$fo_STARTING_FASTA".filtered.bak

tr '-' '_' < "$fo_STARTING_FASTA".filtered.bak | sed 's/_L001/-L001/g' > "$fo_STARTING_FASTA".filtered


echo
echo "De novo chimera detection"

$VSEARCH --threads "$THREADS" \
    --uchime_denovo "$fo_STARTING_FASTA".filtered \
    --sizein \
    --fasta_width "$FASTA_WIDTH" \
    --nonchimeras "$fo_ALL_DENOVO_NONCHIMERAS_FASTA"


echo "Sequences remaining after de novo chimera detection: $(grep -c "^>" "$fo_ALL_DENOVO_NONCHIMERAS_FASTA")"


echo
echo "Reference chimera detection"

$VSEARCH --threads "$THREADS" \
    --uchime_ref "$fo_ALL_DENOVO_NONCHIMERAS_FASTA" \
    --db "$fi_TAXDB_FASTA" \
    --sizein \
    --fasta_width "$FASTA_WIDTH" \
    --nonchimeras "$fo_ALL_NONCHIMERAS_FASTA"

echo "Sequences remaining after reference-based chimera detection: $(grep -c "^>" "$fo_ALL_NONCHIMERAS_FASTA")"




echo
echo "Cluster at "$IDENT" identity, method="$CLUSTER_METHOD", iddef="$IDDEF" and relabel with OTU_n, generate OTU table"

$VSEARCH --threads "$THREADS" \
    --$CLUSTER_METHOD "$fo_ALL_NONCHIMERAS_FASTA" \
    --id "$IDENT" \
    --iddef "$IDDEF" \
    --userfields "id4" \
    --clusterout_id \
    --strand both \
    --sizein \
    --sizeout \
    --fasta_width "$FASTA_WIDTH" \
    --uc "$fo_CLUSTERING_UC" \
    --relabel OTU_ \
    --centroids "$fo_ALL_OTUS_FASTA" \
    --maxaccepts 0 \
    --maxrejects 0 \
    --otutabout "$fo_ALL_OTUTAB_TXT"
    

echo
echo "Number of OTUs: $(grep -c "^>" "$fo_ALL_OTUS_FASTA")"


echo
echo "Taxonomic classification, sintax_cutoff="$SINTAX_CUTOFF""

$VSEARCH --threads "$THREADS" \
    --sintax "$fo_ALL_OTUS_FASTA" \
    --db "$fi_TAXDB_FASTA" \
    --sintax_cutoff "$SINTAX_CUTOFF" \
    --strand both \
    --tabbedout "$fo_ALL_TAXED_OUT"



echo
echo "Filling in taxonomy for OTU-abundance table..." #10 mins for 8000 OTUs previously


sort -V "$fo_ALL_TAXED_OUT" | awk '{print $1,$4}' | sed -E 's/(OTU_[0-9]{1,})/\1#/g' | tr ' ' '#' > "$fo_ALL_TAXED_OUT".tmp

tail -n +2 "$fo_ALL_OTUTAB_TXT" | sort -V | sed -E 's/(OTU_[0-9]{1,})/\1#/g' | tr '\t' ';' > "$fo_ALL_OTUTAB_TXT".tmp

awk -F'#' '{print $1}' "$fo_ALL_TAXED_OUT".tmp > "$fo_ALL_TAXED_OUT".list.tmp
awk -F'#' '{print $1}' "$fo_ALL_OTUTAB_TXT".tmp > "$fo_ALL_OTUTAB_TXT".list.tmp

difftest=$(diff "$fo_ALL_TAXED_OUT".list.tmp "$fo_ALL_OTUTAB_TXT".list.tmp)

if [ "$difftest" == "" ]
then
	join -t'#' -1 1 -2 1 -o 1.1,1.2,1.3,2.2 "$fo_ALL_TAXED_OUT".tmp "$fo_ALL_OTUTAB_TXT".tmp > "$fo_COUNTS_WITH_TAX".nolabel.tmp

else
	echo  "Lists do not match, cannot perform join. Exiting." $COL_RESET
	exit
fi

sed 's/##/#Unassigned#/g' "$fo_COUNTS_WITH_TAX".nolabel.tmp > "$fo_COUNTS_WITH_TAX".tmp

head -n 1 "$fo_ALL_OTUTAB_TXT" | sed 's/#/TAX;/' | sed 's/OTU ID/OTU_ID/g' | tr ';' '\t' > "$fo_COUNTS_WITH_TAX"
awk -F'#' '{OFS="\t"; print $3,$1,$4}' "$fo_COUNTS_WITH_TAX".tmp | sed '/^;/d' | tr ';' '\t' | tr -s '\t' '\t' >> "$fo_COUNTS_WITH_TAX"


echo
echo "Aggregating counts..."

Rscript "$SCRIPTDIR"/aggregate.R "$fo_COUNTS_WITH_TAX" "$fo_AGG_COUNTS"


rm "$OUTDIR"/*.tmp
echo  "Done." $COL_RESET
date
exit


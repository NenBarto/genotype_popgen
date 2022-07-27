
conda activate enviroDNA

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/results/1905CR5"

#amplicons=( "18S" "COI" )
amplicons=( "COI" )
types=( "_all_0.9" "_pests_0.98" )

for amplicon in ${amplicons[@]};do
  for type in ${types[@]}; do
    annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/$amplicon"
    #blastLine="blastn -query $inDir/all.otus_$amplicon$type.fasta -db $annotationDir/"$amplicon"_metazoa -outfmt 6 -out all.otus-"$amplicon$type"_metazoa.tab"
  	blastLine="blastn -word_size 11 -qcov_hsp_perc 50 -query $inDir/all.otus_$amplicon$type.fasta -db $annotationDir/"$amplicon"_metazoa_uniq -outfmt \"6 qseqid sseqid pident length mismatch gapopen qcovs evalue bitscore\" -out all.otus-"$amplicon$type"_metazoa.tab"

    qsub -b y -cwd -j y -N blast$amplicon$type -R y -pe smp 20 -V $blastLine
  done;
done;


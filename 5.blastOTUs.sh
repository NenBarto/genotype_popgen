
conda activate enviroDNA

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/results/1905CR5/split"
outDir="/share/ScratchGeneral/nenbar/projects/Anthony/results/1905CR5/out"
mkdir -P $outDir
#amplicons=( "18S" "COI" )
amplicons=( "COI" )
types=( "_all_0.9" )

for amplicon in ${amplicons[@]};do
  for type in ${types[@]}; do
	for file in `ls $inDir`; do
	#	echo $file
    annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/$amplicon"
    #blastLine="blastn -query $inDir/all.otus_$amplicon$type.fasta -db $annotationDir/"$amplicon"_metazoa -outfmt 6 -out all.otus-"$amplicon$type"_metazoa.tab"
  	blastLine="blastn -word_size 11 -qcov_hsp_perc 50 -query $inDir/$file -db $annotationDir/"$amplicon"_metazoa_uniq -outfmt \"6 qseqid sseqid pident length mismatch gapopen qcovs evalue bitscore\" -out $outDir/all.otus.$file.tab"
#echo $blastLine
  qsub -b y -cwd -j y -N blast$file -R y -pe smp 5 -V $blastLine
	done;  
done;
done;


conda activate enviroDNA


projectname="COI"
scriptDir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/$projectname"
projectDir=$scriptDir
inDirs="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/$projectname/split"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir
mkdir -p $inDirs


###### first fetch the sequences

export PATH=${PATH}:${HOME}/edirect
export NCBI_API_KEY=60b6e6ed96096ce6ea0f8c049b625ffb6a09
#esearch -db nucleotide -query "(COI AND mitochondrial)) AND 1:10000[Sequence Length] AND animals[filter]" | efetch -format fasta > animals.fasta

#esearch -db nucleotide -query "(18S AND (rRNA OR ribosomal))  AND animals[filter]" | efetch -format fasta > animals.fasta
animalLine="esearch -db nucleotide -query \"(COI AND mitochondrial) AND 1:10000[Sequence Length] AND animals[filter]\" | efetch -format fasta > animals.fasta"
qsub -b y -cwd -j y -N animalTest -R y -pe smp 1 -V $animalLine

metazoa_COI <- entrez_search(db="nuccore", term="COI[base] AND mitochondrial[base] AND metazoa[ORGN]", use_history=TRUE)

#split the search term in ten intervals
#each

#for( seq_start in seq(i*100000+1,c(i+1)*100000,250)){
#  cat(seq_start)
#  cat("\n")
#    recs <- entrez_fetch(db="nuccore", web_history=metazoa_COI$web_history,
#                         rettype="fasta", retmax=250, retstart=seq_start)
#    cat(recs, file=paste0("COI_test.fasta"), append=TRUE)
#    cat(seq_start+249, "sequences downloaded\r")
#}


###### WARNING in this case we used instead the premade database from https://www.nature.com/articles/sdata2018156
###### file CO-ARBitrator_joint.fasta was made by merging the first version and the increment
#prepare file for pyfasta
perl -pe 'print $_;$i=0;while(<>){$i++;if($_=~m/>/){s/\n/.$i\n/;};print $_}' <temp.fasta >temp.uniq.fasta

#split fasta files with pyfasta
pyfasta split -n 10 temp.uniq.fasta

mv *0*.fasta $inDirs


for file in `ls $inDirs/*.fasta`;do
  echo $file
  pyfasta split -n 20 $file
done;

#mv $inDirs/CMetaCOXI_Seqs.0{0..9}.fasta ../seqs

for file in `ls $inDirs/*.fasta | grep -P "\\d.fasta"`;do
  sampleName=`basename $file | sed 's/.fasta//'`
  mkdir -p $inDirs/$sampleName
  mv -f $file $inDirs/$sampleName
done

#rm $inDirs/*fasta* 


for inDir in `ls $inDirs `;do
  sampleName=`echo $inDir | sed 's/\///'`
  echo $sampleName
  hammingLine="python $scriptDir/hammingdapt4.py --input_dir $inDirs/$sampleName -f GGWACWGGWTGAACWGTWTAYCCYCC -r TAIACYTCIGGRTGICCRAARAAYCA --min_length=60 --max_length=1000 --distance_limit=15"
  echo $hammingLine
  qsub -b y -cwd -j y -N hd -R y -pe smp 2 -V $hammingLine
done

#hardy: 131330
for file in `ls $inDirs/**/*.trimmed.fasta`;do
  sampleName=`echo $file | sed 's/.trimmed.fasta/.coi.fasta/'`
  echo $sampleName
  cp $file $sampleName
done

mkdir db
mv $inDirs/**/*coi* db

for inDir in `ls $inDirs`;do
  sampleName=`echo $inDir | sed 's/\///'`
  echo $sampleName
  hammingLine="python $scriptDir/hammingdapt4.py --input_dir $inDirs/$sampleName -f AGGGCAAKYCTGGTGCCAGC -r GRCGGTATCTRATCGYCTT --min_length=60 --max_length=1000 --distance_limit=7"
  qsub -b y -cwd -j y -N hd -R y -pe smp 1 -V $hammingLine
done

for file in `ls $inDirs/**/*.trimmed.fasta`;do
  sampleName=`echo $file | sed 's/.trimmed.fasta/.coi.fasta/'`
  echo $sampleName
  mv $file $sampleName
done

rm $scriptDir/split/**/*trimmed*
mkdir $scriptDir/metaCOXI
mv $scriptDir/split/**/*coi* $scriptDir/metaCOXI
cat $scriptDir/metaCOXI/* >metaCOXI.fasta

for inDir in `ls $inDirs`;do
  sampleName=`echo $inDir | sed 's/\///'`
  echo $sampleName
  hammingLine="python $scriptDir/hammingdapt4.py --input_dir $inDirs/$sampleName -f CCAGCASCYGCGGTAATTCC -r ACTTTCGTTCTTGATYRATGA --min_length=60 --max_length=1000 --distance_limit=15"
  qsub -b y -cwd -j y -N hd -R y -pe smp 1 -V $hammingLine
done

#V4: 117334
for file in `ls $inDirs/**/*.trimmed.fasta`;do
  sampleName=`echo $file | sed 's/.trimmed.fasta/.V4.fasta/'`
  echo $sampleName
  mv $file $sampleName
done

rm $scriptDir/split/**/*trimmed*
mkdir V4
mv $scriptDir/split/**/*V4* V4
cat $scriptDir/V4/* >V4.fasta







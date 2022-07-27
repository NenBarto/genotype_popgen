
#originally there are 1276458 sequences

conda init 
conda activate enviroDNA

projectname="18S"
scriptDir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/$projectname"
projectDir=$scriptDir
inDirs="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/$projectname/split"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

dbDir="$projectDir/db/pests"


#split fasta files with pyfasta
pyfasta split -n 10 $dbDir/pests_uniq.fasta

mv $dbDir/*0*.fasta $inDirs



for file in `ls $inDirs/*.fasta | grep -P "\\d.fasta"`;do
  sampleName=`basename $file | sed 's/.fasta//'`
  mkdir -p $inDirs/$sampleName
  mv -f $file $inDirs/$sampleName
done

rm $inDirs/*fasta* 


cutoff=11


for inDir in `ls $inDirs`;do
  sampleName=`echo $inDir | sed 's/\///'`
  echo $sampleName
  hammingLine="python $scriptDir/hammingdapt4.py --input_dir $inDirs/$sampleName -f AGGGCAAKYCTGGTGCCAGC -r GRCGGTATCTRATCGYCTT --min_length=60 --max_length=1000 --distance_limit=$cutoff"
  qsub -b y -cwd -j y -N hd -R y -pe smp 2 -V $hammingLine
done


for file in `ls $inDirs/**/*.trimmed.fasta`;do
  sampleName=`echo $file | sed 's/.trimmed.fasta/.pests.fasta/'`
  echo $sampleName
  mv -f $file $sampleName
done

rm $scriptDir/split/**/*trimmed*
mkdir pests.$cutoff
mv -f $scriptDir/split/**/*pests.fasta pests.$cutoff
cat $scriptDir/pests.$cutoff/* >pests.$cutoff.fasta

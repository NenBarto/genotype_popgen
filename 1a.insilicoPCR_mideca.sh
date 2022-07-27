
conda init 
conda activate enviroDNA

projectname="16S"
primers="mideca"
scriptDir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/16S"
inDirs="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/16S/split"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results/"$primers
mkdir -p $resultsDir


for inDir in `ls $inDirs | grep 16S`;do
  sampleName=`echo $inDir | sed 's/\///'`
  hammingLine="python $scriptDir/hammingdapt4.py --input_dir $inDirs/$sampleName -f GGACGATAAGACCCTATAAA -r ACGCTGTTATCCCTAAAGT --min_length=70 --max_length=600 --distance_limit=5 --include_primers"
  qsub -b y -cwd -j y -N hd -R y -pe smp 1 -V $hammingLine
done

moveLine="mv $inDirs/**/*trimmed* $resultsDir"
qsub -b y -cwd -j y -hold_jid hd -N mv -R y -pe smp 1 -V $moveLine


#this results in 32406 untrimmed and
#125875 trimmed
#now where are the rest
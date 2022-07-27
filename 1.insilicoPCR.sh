
conda init 
conda activate enviroDNA

projectname="16S"
scriptDir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/16S"
inDirs="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/16S/split"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir


for inDir in `ls $inDirs | grep 16S`;do
  sampleName=`echo $inDir | sed 's/\///'`
  hammingLine="python $scriptDir/hammingdapt4.py --input_dir $inDirs/$sampleName -f RGACGAGAAGACCCTATARA -r ACGCTGTTATCCCTAARGTA --min_length=70 --max_length=600 --distance_limit=5 --include_primers"
  qsub -b y -cwd -j y -N hd$sampleName -R y -pe smp 1 -V $hammingLine
done

#this results in 32406 untrimmed and
#125875 trimmed
#now where are the rest
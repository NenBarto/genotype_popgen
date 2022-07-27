#conda init 
#conda activate gatk


#to do
#GATK
#R package

projectname="EBB"

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/EBB"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/EBB"
index=$annotationDir/"EBB.scaffolds"
logDir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/EBB/logs"
mkdir -p $logDir




trimDir="$resultsDir/"$projectname".trimgalore"
mkdir -p $trimDir
bwaDir="$resultsDir/"$projectname".bwa"
mkdir -p $bwaDir


#load files
files=`ls $inDir/*R1.fastq.gz`
files=( $files )

for file in ${files[@]}; do 
  echo $file; 
  #file="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq/7111_S29_L001_R1_001.fastq.gz"

  sampleName=`basename $file | sed s/_R1.fastq.gz//`
  echo $sampleName

  file1=$inDir/$sampleName"_R1.fastq.gz"
  file2=$inDir/$sampleName"_R2.fastq.gz"
  fileTrim1gz=$trimDir/$sampleName"_R1_001_val_1.fq.gz"
  fileTrim2gz=$trimDir/$sampleName"_R2_001_val_2.fq.gz"
  bamFile="$bwaDir/$sampleName.bam"
  sortedBamFile="$bwaDir/$sampleName.sorted"
  samFile="$bwaDir/$sampleName.sam"

  trim_line="trim_galore $file1 $file2 --fastqc --paired --retain_unpaired --length 16 -o $trimDir"

  #set up header for read groups, important or breaks GATK
  HEADER=`printf @RG%sID:%s%sSM:%s%sPL:ILLUMINA '\\t' $sampleName '\\t' $sampleName '\\t'`
  mapping_line="bwa mem -R '"${HEADER}"' $index $fileTrim1gz $fileTrim2gz -o $samFile" 
  sam2bam_line="samtools view -S -b $samFile -o $bwaDir/$sampleName.bam"
  sort_line="samtools sort $bwaDir/$sampleName.bam $sortedBamFile"
  index_line="samtools index $bwaDir/$sampleName.sorted.bam"
  clean_line="rm -f $bamFile & rm -f $samFile"

  #qsub -b y -hold_jid indexGenome -wd $logDir -j y -N trim$sampleName -R y -pe smp 1 -V $trim_line
  qsub -b y -hold_jid trim_line$sampleName -wd $logDir -j y -N map$sampleName -R y -pe smp 5 -V $mapping_line
  qsub -b y -hold_jid map$sampleName -wd $logDir -j y -N sam2bam$sampleName -R y -pe smp 12 -V $sam2bam_line
  qsub -b y -hold_jid sam2bam$sampleName -wd $logDir -j y -N sort$sampleName -R y -pe smp 5 -V $sort_line
  qsub -b y -hold_jid sort$sampleName -wd $logDir -j y -N index$sampleName -R y -pe smp 1 -V $index_line
  qsub -b y -hold_jid index$sampleName -wd $logDir -j y -N clean$sampleName -R y -pe smp 1 -V $clean_line

done
#map files


#IGV
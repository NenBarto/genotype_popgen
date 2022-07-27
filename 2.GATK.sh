conda init 
conda activate enviroDNA
module load centos6.10/gi/java/jdk1.8.0_25
#source activate gatk - could be potentially used instead of the containers

#to do
#GATK
#mark duplicates DONE
##Variant Quality Score Recalibration (VQSR) DONE
#integrate GVCF files 
#genotype GVCFs
#R package

projectname="EBB"

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/$projectname/"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/EBB"
index=$annotationDir"/EBB.scaffolds.fa"

logDir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/EBB/logs"
mkdir -p $logDir

#trimming
trimDir="$resultsDir/"$projectname".trimgalore"
mkdir -p $trimDir
bwaDir="$resultsDir/"$projectname".bwa"
mkdir -p $bwaDir
mrkdupDir="$resultsDir/"$projectname".md-alignments"
mkdir -p $mrkdupDir

cutadapt_minlen=30

tempGVCFdir="$resultsDir/"$projectname".temp-gvcf"
mkdir -p $tempGVCFdir



FWD="CTGGTGATTCCCTGTGACG"
REV="AATGTCTCCACCGCGTTC"
REVrc="GAACGCGGTGGAGACATT"
e1=0.2 #accept 20% error rate in primer sequence
e2=0.2
cutadapt_minlen=30
REFERENCE="EBB.scaffolds.fa"
THREADS=5
GATKimage="/share/ScratchGeneral/nenbar/local/images/gatk_latest.sif"

#load files
files=`ls $inDir/*R1.fastq.gz | grep -v F1.F1`
files=( $files )

for file in ${files[@]}; do 
  echo $file; 
  #file="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq/7111_S29_L001_R1_001.fastq.gz"

  sampleName=`basename $file | sed s/_R1.fastq.gz//`
  echo $sampleName
  bamFile="$bwaDir/$sampleName.bam"

  #perform mapping from 1.map.sh 
  ####### notice the change of the name and fix it for the sorted bam
  mappedBamFile="$sampleName.sorted.bam"
  MDfile1=$sampleName.md.bam
  MDfile2=$sampleName.md.fx.bam
  MDfile2temp=$sampleName.md.fx.bam.temp
  VQSRfile1=$sampleName.g.vcf.gz
  reheaderCommand=$logDir/$sampleName".reheader.txt"

  mem="4G"
#  nCores=5

  #md_line0="samtools view -H $bwaDir/$mappedBamFile | grep -v -f $annotationDir/backup/nullContigs.txt > $bwaDir/template2$sampleName.sam; samtools reheader $bwaDir/template2$sampleName.sam $bwaDir/$mappedBamFile > $bwaDir/$mappedBamFile.temp; cp $bwaDir/$mappedBamFile.temp $bwaDir/$mappedBamFile"

  md_line1="singularity exec --bind $bwaDir:/bwaDir,$mrkdupDir:/mrkdupDir $GATKimage gatk --java-options "-Xmx$mem" MarkDuplicates -I /bwaDir/$mappedBamFile -O /mrkdupDir/$MDfile1 -M /mrkdupDir/$sampleName.txt -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000"
  md_line2="singularity exec --bind $bwaDir:/bwaDir,$mrkdupDir:/mrkdupDir,$annotationDir:/annotationDir $GATKimage gatk --java-options "-Xmx$mem" SetNmMdAndUqTags --INPUT /mrkdupDir/$MDfile1 --OUTPUT /mrkdupDir/$MDfile2 --REFERENCE_SEQUENCE /annotationDir/$REFERENCE"
  md_line3="singularity exec --bind $bwaDir:/bwaDir,$mrkdupDir:/mrkdupDir $GATKimage gatk --java-options "-Xmx$mem" BuildBamIndex --INPUT /mrkdupDir/$MDfile2"
  VQSR_line1="singularity exec --bind $mrkdupDir:/mrkdupDir,$tempGVCFdir:/tempGVCFdir,$annotationDir:/annotationDir $GATKimage gatk --java-options \"-Xms20G -Xmx20G -XX:ParallelGCThreads=2\" HaplotypeCaller -R /annotationDir/$REFERENCE -I /mrkdupDir/$MDfile2 -O /tempGVCFdir/$VQSRfile1 -ERC GVCF --native-pair-hmm-threads 8"
  rm_line="rm $mrkdupDir/$MDfile1 && rm $mrkdupDir/$MDfile1.bai && rm $mrkdupDir/$MDfile1.sbi &&"



  echo $md_line0 > $reheaderCommand
  chmod 775 $reheaderCommand

  #qsub -b y -hold_jid clean$sampleName -wd $logDir -j y -N md0$sampleName -R y -pe smp $THREADS -V $reheaderCommand
  #qsub -b y -hold_jid clean$sampleName -wd $logDir -j y -N md1$sampleName -R y -pe smp $THREADS -V $md_line1
  #qsub -b y -hold_jid md1$sampleName -wd $logDir -j y -N md2$sampleName -R y -pe smp $THREADS -V $md_line2
  #qsub -b y -hold_jid md2$sampleName -wd $logDir -j y -N md3$sampleName -R y -pe smp $THREADS -V $md_line3
  #for i in {1..10};do
  #  VQSR_line1="singularity exec --bind $mrkdupDir:/mrkdupDir,$tempGVCFdir:/tempGVCFdir,$annotationDir:/annotationDir $GATKimage gatk --java-options \"-Xms20G -Xmx20G -XX:ParallelGCThreads=2\" HaplotypeCaller -R /annotationDir/$REFERENCE -L /annotationDir/EBB_chrs_$i.bed -I /mrkdupDir/$MDfile2 -O /tempGVCFdir/$VQSRfile1.$i -ERC GVCF --native-pair-hmm-threads 8"    
  #  qsub -b y -hold_jid md3$sampleName -wd $logDir -j y -N VQSR$sampleName -R y -pe smp $THREADS -V $VQSR_line1;
  #done


  vcfFiles=`ls $tempGVCFdir | grep $sampleName | grep -v idx | grep -v tbi`
  vcfFiles=( $vcfFiles )
  printf -v var "I=%s " "${vcfFiles[@]}"
  sortVCFline="java -jar ~/local/bin/picard.jar SortVcf $var O=$sampleName.sorted.vcf.gz"
  qsub -b y -hold_jid VQSR$sampleName -wd $tempGVCFdir -j y -N sortVCF -R y -pe smp 9 -V $sortVCFline
  #qsub -b y -hold_jid VQSR$sampleName -wd $logDir -j y -N clean2$sampleName -R y -pe smp 1 -V $rm_line

done;

#eliminate the null contigs
#split haplotype calling into 10

library(rtracklayer)
data=read.table("../../annotation/EBB/EBB.scaffolds.fa.fai")
splitList<-cut(1:dim(data)[1],100)
gr<-GRanges(seqnames=data$V1,IRanges(1,data$V2))
#in the future do this to allow for equal size of vcf files and fewer failed processes
gr<-gr[sample(1:length(gr))]
grL<-split(gr,splitList)

grLE<-endoapply(grL,sort)
for(i in 1:length(grLE)){
  export.bed(grLE[[i]],paste0("../../annotation/EBB/EBB_chrs_",i,".bed"))
}

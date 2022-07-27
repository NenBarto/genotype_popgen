
conda init 
conda activate gatk

projectname="210813_miseq_AR"

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation"
index=$annotationDir/"mhc282"
logDir="/share/ScratchGeneral/nenbar/projects/Anthony/logs"
mkdir -p $logDir

cutadapt_minlen=30

#trimming
revcompDir="$resultsDir/"$projectname".revcomp"
mkdir -p $revcompDir
filterDir="$resultsDir/"$projectname".filter"
mkdir -p $filterDir
mergedDir="$resultsDir/"$projectname".merged"
mkdir -p $mergedDir
cutadaptDir="$resultsDir/"$projectname".cutadapt"
mkdir -p $cutadaptDir
dedupDir="$resultsDir/"$projectname".dedup"
mkdir -p $dedupDir
bwaDir="$resultsDir/"$projectname".bwa"
mkdir -p $bwaDir
haplotypesDir="$resultsDir/"$projectname".haplotypes"
mkdir -p $haplotypesDir
tempRecabBamDir="$resultsDir/"$projectname".tempRecabBam"
mkdir -p $tempRecabBamDir
recabBamDir="$resultsDir/"$projectname".recabBam"
mkdir -p $recabBamDir

FWD="CTGGTGATTCCCTGTGACG"
REV="AATGTCTCCACCGCGTTC"
REVrc="GAACGCGGTGGAGACATT"
e1=0.2 #accept 20% error rate in primer sequence
e2=0.2
cutadapt_minlen=30
THREADS=5
MINOVLEN=10 #overlap length; default 10
MAXDIFFPCT=100 #overlap cannot have % mismatches higher than this threshold; default 100
MAXEE=2 #Elbrecht 2016 uses maxee=1
MAXEE2=2
MAXDIFFS=10
minuniquesize=10
#load files
files=`ls $inDir/*R1_001.fastq.gz`
files=( $files )


#DBI import
subset="all"
rm $sampleDir/gvcfs-for-db-import.sample_map
for VCF in $tempGVCFdir/*.gz; do
    filename=$(basename -- "$VCF");
    name="${filename%%.*}";
    file=`basename $VCF`
    VCF="/tempGVCFdir/"$file""
    echo -e "$name\t$VCF" >> $sampleDir/gvcfs-for-db-import.sample_map;
done

#make a database of SNPs  
#FIX TEMP!!!  
dbiImport_line="singularity exec --bind $tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir $GATKimage gatk --java-options "-Xmx58g" \
    GenomicsDBImport \
    --genomicsdb-workspace-path subset_$subset \
    --tmp-dir tmp \
    --batch-size 20 \
    --sample-name-map /sampleDir/gvcfs-for-db-import.sample_map \
    --L mhc282:1-282"
$dbiImport_line

# call cohort SNPs
genotype_line="singularity exec --bind $vcfDir:/vcfDir,$annotationDir:/annotationDir,$tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir $GATKimage gatk --java-options "-Xmx58g" \
    GenotypeGVCFs \
    -R /annotationDir/$REFERENCE \
    -V gendb://subset_$subset \
    -O /vcfDir/subset_$subset.vcf.gz"
$genotype_line

#perform vcftools analysis
vcftools --gzvcf $vcfDir/subset_$subset.vcf.gz --indv-freq-burden --out $resultDir/subset_$subset
vcftools --gzvcf $vcfDir/subset_$subset.vcf.gz --freq --out $resultDir/subset_$subset

#vcftools
#vcftools \
#    --gzvcf $vcfDir/subsample_711.vcf.gz \
#    --minDP 30 \
#    --minQ 30 \
#    --minGQ 30 \
#    --remove-indels --out $resultDir/subsample_711


for file in ${files[@]}; do 
  echo $file; 
  #file="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq/7111_S29_L001_R1_001.fastq.gz"

  sampleName=`basename $file | sed s/_L001_R1_001.fastq.gz//`
  echo $sampleName

  INPUT=$sampleName.sorted.bam
  TEMP_BQSR=$sampleName.bqsr-recal
  OUTBQSR=$sampleName.BQSR.bam
  
  sortedBamFile="$bwaDir/$sampleName.sorted.bam"

  baseRecab1_line="gatk --java-options "-Xmx4G" BaseRecalibrator \
    -I /bwaDir/$INPUT \
    -R/annotationDir/$REFERENCE \
    --known-sites /genotype/$SITES \
    -O /tempRecabBam/TEMP_BQSR"
  baseRecab2_line="gatk --java-options "-Xmx4G" ApplyBQSR \
    -R $REFERENCE \
    -I /bwaDir/$INPUT \
    --bqsr-recal-file /tempRecabBam/TEMP_BQSR \
    -O /recabBam/$OUTBQSR"
  qsub -b y -hold_jid clean$sampleName -wd $logDir -j y -N baseRecab1 -R y -pe smp $THREADS -V $baseRecab1_line
  qsub -b y -hold_jid baseRecab1 -wd $logDir -j y -N baseRecab2 -R y -pe smp $THREADS -V $baseRecab2_line

done;
  



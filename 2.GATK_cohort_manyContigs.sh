
conda init 
conda activate gatk

projectname="EBB"

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/$projectname/"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/$projectname"
index=$annotationDir"/EBB.scaffolds.fa"

logDir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/$projectname/logs"
mkdir -p $logDir

GATKimage="/share/ScratchGeneral/nenbar/local/images/gatk_latest.sif"

cutadapt_minlen=30

sampleDir="$resultsDir/"$projectname".samples"
mkdir -p $sampleDir
tempGVCFdir="$resultsDir/"$projectname".temp-gvcf"
mkdir -p $tempGVCFdir
tmpDir="$resultsDir/"$projectname".temp"
mkdir -p $tmpDir
vcfDir="$resultsDir/"$projectname".vcf"
mkdir -p $vcfDir
resultDir="$resultsDir/"$projectname".result"
mkdir -p $resultDir
filteredDir="$resultsDir/"$projectname".vcffiltered"
mkdir -p $filteredDir
nCores=10

REFERENCE="EBB.scaffolds.fa"

#load files
files=`ls $inDir/*R1.fastq.gz`
files=( $files )


#DBI import

rm $sampleDir/gvcfs-for-db-import.sample_map
for VCF in $tempGVCFdir/*sorted*.gz; do
    filename=$(basename -- "$VCF");
    name="${filename%%.sorted.*}";
    file=`basename $VCF`
    VCF="/tempGVCFdir/"$file""
    echo -e "$name\t$VCF" >> $sampleDir/gvcfs-for-db-import.sample_map;
done

#make a database of SNPs but in mind with large number of contigs use
#--genomicsdb-shared-posixfs-optimizations
#


for i in {1..100};do

    subset="chr$i"
    dbiImport_line="singularity exec --bind $tmpDir:/tmpDir,$tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir,$annotationDir:/annotationDir $GATKimage gatk --java-options "-Xmx60g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path $subset \
        --tmp-dir /tmpDir \
        --sample-name-map /sampleDir/gvcfs-for-db-import.sample_map \
        --merge-contigs-into-num-partitions 10 \
        --genomicsdb-shared-posixfs-optimizations true \
        --L /annotationDir/EBB_chrs_$i.bed"

    qsub -b y -hold_jid sortVCF -wd $logDir -j y -N dbiImport$subset -R y -pe smp $nCores -V $dbiImport_line

    # call cohort SNPs
    genotype_line="singularity exec --bind $vcfDir:/vcfDir,$annotationDir:/annotationDir,$tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir $GATKimage gatk --java-options "-Xmx60g" \
        GenotypeGVCFs \
        -R /annotationDir/$REFERENCE \
        -V gendb://$subset \
        -O /vcfDir/$subset.vcf.gz"
    qsub -b y -hold_jid dbiImport$subset -wd $logDir -j y -N genotype -R y -pe smp $nCores -V $genotype_line
done;

#the following filtering is done based on the ovarflow workflow
#https://ovarflow.readthedocs.io/en/latest/TheWorkflow/Overview.html?highlight=filtering#hard-filtering

#library(rtracklayer)
#data=read.table("../../annotation/EBB/EBB.scaffolds.fa.fai")
#splitList<-cut(1:dim(data)[1],100)
#gr<-GRanges(seqnames=data$V1,IRanges(1,data$V2))
##in the future do this to allow for equal size of vcf files and fewer failed processes
#gr<-gr[sample(1:length(gr))]
#grL<-split(gr,splitList)
#
#grLE<-endoapply(grL,sort)
#for(i in 1:length(grLE)){
#  export.bed(grLE[[i]],paste0("EBB_chrs_",i,".bed"))
#}



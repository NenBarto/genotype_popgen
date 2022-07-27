
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
vcfDir="$resultsDir/"$projectname".vcf"
mkdir -p $vcfDir
resultDir="$resultsDir/"$projectname".result"
mkdir -p $resultDir


REFERENCE="EBB.scaffolds.fa"

#load files
files=`ls $inDir/*R1.fastq.gz`
files=( $files )


#DBI import
subset="mitochondrial"
rm $sampleDir/gvcfs-for-db-import.sample_map
for VCF in $tempGVCFdir/*sorted*.gz; do
    filename=$(basename -- "$VCF");
    name="${filename%%.sorted.*}";
    file=`basename $VCF`
    VCF="/tempGVCFdir/"$file""
    echo -e "$name\t$VCF" >> $sampleDir/gvcfs-for-db-import.sample_map;
done

#make a database of SNPs   
dbiImport_line="singularity exec --bind $tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir $GATKimage gatk --java-options "-Xmx58g" \
    GenomicsDBImport \
    --genomicsdb-workspace-path $subset \
    --tmp-dir tmp \
    --batch-size 20 \
    --sample-name-map /sampleDir/gvcfs-for-db-import.sample_map \
    --L scaffold61811"

  qsub -b y -hold_jid sortVCF -wd $logDir -j y -N dbiImport -R y -pe smp 20 -V $dbiImport_line

# call cohort SNPs
genotype_line="singularity exec --bind $vcfDir:/vcfDir,$annotationDir:/annotationDir,$tempGVCFdir:/tempGVCFdir,$sampleDir:/sampleDir $GATKimage gatk --java-options "-Xmx58g" \
    GenotypeGVCFs \
    -R /annotationDir/$REFERENCE \
    -V gendb://initial \
    -O /vcfDir/$subset.vcf.gz"
qsub -b y -hold_jid dbiImport -wd $logDir -j y -N genotype -R y -pe smp 20 -V $genotype_line

#perform vcftools analysis
#vcftools --gzvcf $vcfDir/$subset.vcf.gz --indv-freq-burden --out $resultDir/$subset
#vcftools --gzvcf $vcfDir/$subset.vcf.gz --freq --out $resultDir/$subset
#
##vcftools
#vcftools \
#    --gzvcf $vcfDir/$subset.vcf.gz \
#    --minDP 30 \
#    --minQ 30 \
#    --minGQ 30 \
#    --remove-indels --out $resultDir/$subset.filtered



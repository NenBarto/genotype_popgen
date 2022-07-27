
conda init 
conda activate gatk

projectname="EBB"


#### DIRECTORIES
inDir="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/$projectname/"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/$projectname"
index=$annotationDir"/EBB.scaffolds.fa"

logDir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/$projectname/logs"
mkdir -p $logDir

vcfDir="$resultsDir/"$projectname".vcf"
mkdir -p $vcfDir
filteredDir="$resultsDir/"$projectname".vcffiltered"
mkdir -p $filteredDir

#### FILES

GATKimage="/share/ScratchGeneral/nenbar/local/images/gatk_latest.sif"

sortedFile=$filteredDir/$projectname.sorted.vcf.gz
sortedFileSNPs=$filteredDir/$projectname.SNPs.sorted.vcf.gz
sortedFileINDs=$filteredDir/$projectname.INDs.sorted.vcf.gz

#### PARAMETERS

cutadapt_minlen=30
nCores=10


#### ANALYSIS


#combine and sort all the vcfs
vcfFiles=`ls $vcfDir | grep chr | grep -v tbi`
vcfFiles=( $vcfFiles )
printf -v var "I=$vcfDir/%s " "${vcfFiles[@]}"
sortVCFline="java -jar ~/local/bin/picard.jar SortVcf $var O=$sortedFile"
qsub -b y -hold_jid genotype -wd $logDir -j y -N sortVCF -R y -pe smp $nCores -V $sortVCFline

#separate into SNPs and indels
separateSNPline="singularity exec --bind $vcfDir:/vcfDir,$filteredDir:/filteredDir $GATKimage gatk --java-options "-Xmx60g" \
    SelectVariants -V /filteredDir/$projectname.sorted.vcf.gz -select-type SNP -O /filteredDir/$projectname.SNPs.sorted.vcf.gz"
qsub -b y -hold_jid sortVCF -wd $logDir -j y -N getSNPs -R y -pe smp $nCores -V $separateSNPline

separateINDELsline="singularity exec --bind $vcfDir:/vcfDir,$filteredDir:/filteredDir $GATKimage gatk --java-options "-Xmx60g" \
    SelectVariants -V /filteredDir/$projectname.sorted.vcf.gz -select-type INDEL -select-type MIXED -O /filteredDir/$projectname.INDs.sorted.vcf.gz"
qsub -b y -hold_jid sortVCF -wd $logDir -j y -N getINDs -R y -pe smp $nCores -V $separateINDELsline

#filter SNPs

snp_line="singularity exec --bind $vcfDir:/vcfDir,$filteredDir:/filteredDir $GATKimage gatk --java-options "-Xmx60g" \
        VariantFiltration -V /filteredDir/$projectname.SNPs.sorted.vcf.gz \
      -filter 'QD <2.0' --filter-name 'QD2' \
      -filter 'QUAL < 30.0' --filter-name 'QUAL30' \
      -filter 'SOR > 3.0' --filter-name 'SOR3' \
      -filter 'FS > 60.0' --filter-name 'FS60' \
      -filter 'MQ < 40.0' --filter-name 'MQ40' \
      -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' \
      -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' \
      -O /filteredDir/$projectname.SNPs.filt.vcf.gz"
qsub -b y -hold_jid getSNPs -wd $logDir -j y -N filterSNP -R y -pe smp $nCores -V $snp_line



#filter indels
indel_line="singularity exec --bind $vcfDir:/vcfDir,$filteredDir:/filteredDir $GATKimage gatk --java-options "-Xmx18g" \
    VariantFiltration -V /filteredDir/$projectname.INDs.sorted.vcf.gz \
  -filter 'QD < 2.0' --filter-name 'QD2' \
  -filter 'QUAL < 30.0' --filter-name 'QUAL30' \
  -filter 'FS > 200.0' --filter-name 'FS200' \
  -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' \
  -O /filteredDir/$projectname.INDs.filt.vcf.gz"
qsub -b y -hold_jid getINDs -wd $logDir -j y -N filterSNP -R y -pe smp $nCores -V $indel_line

#combine and sort
combine_line="java -jar ~/local/bin/picard.jar SortVcf I=$filteredDir/$projectname.SNPs.filt.vcf.gz I=$filteredDir/$projectname.INDs.filt.vcf.gz  \
  O=$filteredDir/$projectname.sorted_variants.vcf.gz"
qsub -b y -hold_jid filterSNP -wd $logDir -j y -N combineAll -R y -pe smp $nCores -V $combine_line

#hard filtering
hard_filter="singularity exec --bind $vcfDir:/vcfDir,$filteredDir:/filteredDir $GATKimage gatk --java-options "-Xmx18g" \
    SelectVariants -V /filteredDir/$projectname.sorted_variants.vcf.gz \
  -O /filteredDir/$projectname.filtered_variants.vcf.gz --exclude-filtered true"
qsub -b y -hold_jid getINDs -wd $logDir -j y -N filterSNP -R y -pe smp $nCores -V $hard_filter



  

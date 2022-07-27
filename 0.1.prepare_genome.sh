
conda init 
conda activate gatk

projectname="EBB"
annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/EBB/"
REFERENCE="EBB.scaffolds.fa"
#run index of the genome
GATKimage="/share/ScratchGeneral/nenbar/local/images/gatk_latest.sif"
singularity exec --bind $annotationDir:/annotationDir $GATKimage gatk CreateSequenceDictionary -R /annotationDir/$REFERENCE

indexLine="bwa index $annotationDir/EBB.scaffolds.fa"
qsub -b y -cwd -j y -N indexGenome -R y -pe smp 20 -V $indexLine

conda activate R4.0
indir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/COI"

outdir="/share/ScratchGeneral/nenbar/projects/Anthony/scripts/COI/ncbi"

for i in {0..33};do
 echo $i
 nLines=`more $outdir/COI_$i.fasta | grep ">" | wc -l`
 if [ $nLines -lt 100000 ];then
  qsub -b y -cwd -j y -N dwnld$i -R y -pe smp 1 -V "R --vanilla <$indir/0.download.R --param=$i --startnew=$nLines"
 fi
 sleep .27
done;

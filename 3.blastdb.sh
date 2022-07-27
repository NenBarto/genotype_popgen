
conda activate enviroDNA

export BLASTDB_LMDB_MAP_SIZE=100000000

inFa="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/16S_metazoa.fasta.gz"
inTax="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/taxdb.btd"
out="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/16S_metazoa"
commandLine="gunzip -c $inFa | makeblastdb -in - -parse_seqids -blastdb_version 5 -taxid_map $inTax -title 16S_metazoa_261121 -dbtype nucl -out $out"
qsub -b y -cwd -j y -N blast -R y -pe smp 40 -V $commandLine


inFa="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/nt.gz"
inTax="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/taxdb.btd"
out="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/nucleotide"
commandLine="gunzip -c $inFa | makeblastdb -in - -parse_seqids -blastdb_version 5 -taxid_map $inTax -title nucleotide_261121 -dbtype nucl -out $out"
qsub -b y -cwd -j y -N blast -R y -pe smp 40 -V $commandLine

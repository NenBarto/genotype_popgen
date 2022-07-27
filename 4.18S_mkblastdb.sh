conda activate enviroDNA

blastLine="blastn -query all.otus.fasta -db ../../annotation/blast/16S_metazoa -outfmt 6 -out all.otus-16S_metazoa.tab"

export BLASTDB_LMDB_MAP_SIZE=100000000
inFa="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/18S/18S_metazoa.fasta.gz"
inTax="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/18S/taxdb.btd"
out="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/18S/18S_metazoa"
commandLine="gunzip -c $inFa | makeblastdb -in - -parse_seqids -blastdb_version 5 -taxid_map $inTax -title 18S_metazoa_040221 -dbtype nucl -out $out"
qsub -b y -cwd -j y -N blast -R y -pe smp 20 -V $commandLine


inFa="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/nt.gz"
inTax="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/taxdb.btd"
out="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/nucleotide"
commandLine="gunzip -c $inFa | makeblastdb -in - -parse_seqids -blastdb_version 5 -taxid_map $inTax -title nucleotide_261121 -dbtype nucl -out $out"
qsub -b y -cwd -j y -N blast -R y -pe smp 60 -V $commandLine

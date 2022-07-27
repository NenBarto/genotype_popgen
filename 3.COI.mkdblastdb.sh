conda activate enviroDNA

blastLine="blastn -query all.otus.fasta -db ../../annotation/blast/16S_metazoa -outfmt 6 -out all.otus-16S_metazoa.tab"

export BLASTDB_LMDB_MAP_SIZE=100000000
inFa="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/COI/COI_metazoa.short.fasta.gz"
inTax="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/COI/taxdb.btd"
out="/share/ScratchGeneral/nenbar/projects/Anthony/annotation/blast/COI/COI_metazoa_uniq"
commandLine="gunzip -c $inFa | makeblastdb -in - -parse_seqids -blastdb_version 5 -taxid_map $inTax -title COI_metazoa_070221 -dbtype nucl -out $out"
qsub -b y -cwd -j y -N blast -R y -pe smp 20 -V $commandLine



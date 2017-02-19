test $# -lt 1 && {
   echo "usage: $0 contigs.fa"
   exit
}

contigs=$1
# annotate (genes) for contigs ../../longer_contigs/contigs_longer_than_1500bp.fa with prokka
echo "annotate (genes) for contigs $contigs with prokka"

bn=$(basename "$contigs")
dna=${bn%.*}
#"${filename%.*}"
echo $dna

prokka $contigs . --force ON --locustag $dna --genus $dna --species $dna --strain $dna --metagenome ON --cpus 7 
# --cpus: # of cpus.  Defaults to 8. 


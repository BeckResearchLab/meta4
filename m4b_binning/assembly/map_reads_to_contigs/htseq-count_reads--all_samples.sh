echo "use $@: use htseq-count on a all the gffs to get tables of gene counts for each sample"

# safety first
#set -eu 

# commad-line args:
test $# -lt 2 && {
   echo "usage: $0 bam_folder gff_path"
   #echo "e.g. $0 /work/m4b_binning/assembly/data/HOW10/8777.1.112183.AATAGG.fastq.gz /work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa"
   exit
}

alignment_folder=$1
gff_path=$2

bam_files=`ls $1/*/*sorted.bam`
echo "bam files:"
# print the names of the identified bam files
for b in $bam_files
do 
	echo $b 
done

# Function to run htseq-count on one sample: 
# ls /work/m4b_binning/assembly/data/*OW*/*.fastq.gz | parallel --jobs 20 ./map_reads.sh {} /work/m4b_binning/assembly /longer_contigs/contigs_longer_than_1500bp.fa

ls $1/*/*sorted.bam | parallel --jobs 30 ./htseq-count_reads.sh  {} $gff_path 


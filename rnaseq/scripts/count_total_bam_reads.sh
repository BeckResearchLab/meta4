# safety first
# consider -eux for more detailed error reporting. 
set -eu 

echo "Script $0: count total # of reads for each .bam in a set of .bam files."
echo "Warning: $0 over-reports by double-counting reads that map multiple places!!" 

# samtools view -c /work/m4b_binning/analysis/rna-seq/map_to_contigs_longer_than_1500bp/8904.6.115262.TCGGCA/8904.6.115262.TCGGCA.fastq.gz.sorted.bam

# Make sure the script was provided exactly 2 arguments. 
test ! $# -eq 2 && {
   echo "expected 2 arguments; watch out for shell expansion"
   echo "usage: $0 bam_path out_tsv_name"
   echo "example: $0 \"../*/*.sorted.bam\" sample_read_counts.tsv"
   exit
}

# sample bam:
# /work/m4b_binning/analysis/rna-seq/map_to_contigs_longer_than_1500bp/8904.6.115262.TCGGCA/8904.6.115262.TCGGCA.fastq.gz.sorted.bam

BAMS=$1
RESULTS=$2
rm -f $RESULTS

for b in `ls $BAMS`
do
	bam=`basename $b`
	# put the filename in without a new line:
	echo -n "$bam	" >> $RESULTS 
	
	# append the number of @ lines in the file:
	count_command="samtools view -c $b >> $RESULTS "
	echo "count command: $count_command"
	eval $count_command
	
done

echo "done counting all reads in each file"

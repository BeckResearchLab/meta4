# safety first
# consider -eux for more detailed error reporting.
set -eu

echo "Script $0: go into a bunch of alignment folders and make an idxstat file for each.  Index first."

# http://www.htslib.org/doc/samtools.html
# The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.

# Make sure the script was provided exactly 1 argument.
test ! $# -eq 1 && {
   echo "expected 1 arguments; watch out for shell expansion"
   echo "usage: $0 bam_path"
   echo "example: $0 \"../*/*.sorted.bam\""
   exit
}

BAMS=$1

for b in `ls $BAMS`
do
	cd $(dirname "${b}")
	BAM=`basename $b`
	echo bam: $BAM
	OUT_FILE=$BAM.idxstat
	# put the filename in without a new line:
	echo  "contig	bp	mapped	unmapped" > $OUT_FILE

	# append the number of @ lines in the file:
	index_command="samtools index *.bam"
	count_command="samtools idxstats $BAM >> $OUT_FILE "
	echo "index command: $index_command"
	echo "count command: $count_command"
	eval $index_command
	eval $count_command
	cd -  # go back to base dir

done

echo "done tabulating reads assigned to each bam file"

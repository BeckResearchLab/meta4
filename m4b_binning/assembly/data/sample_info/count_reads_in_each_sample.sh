# safety first
set -eu
#set -eux

echo "Script $0: count reads using grep '^@' | wc -l" 

test ! $# -eq 2 && {
#test $# -lt 2 && {
   echo "expected 2 arguments; watch out for shell expansion"
   echo "usage: $0 fastq_path out_tsv_name"
   echo "example: $0 \"../*/*.fastq.gz\" sample_read_counts.tsv"
   exit
}

# sample fastq:
# /work/m4b_binning/assembly/data/HOW11/8777.2.112196.AGATAG.fastq.gz

FASTQS=$1
RESULTS=$2
#RESULTS='sample_read_counts--tmp.tsv'
rm -f $RESULTS

#for fq in `ls ../*/*.fastq.gz`
for fq in `ls $FASTQS`
do
	fastq=`basename $fq`
	#echo $fq
	# put the filename in without a new line:
	echo -n "$fastq	" >> $RESULTS 
	
	# append the number of @ lines in the file:
	count_command="zcat $fq | grep '^@' | wc -l >> $RESULTS "
	echo "count command: $count_command"
	eval $count_command
	
done

echo "done counting all reads in each fastq.gz file"

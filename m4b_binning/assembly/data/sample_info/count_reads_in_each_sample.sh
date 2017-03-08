# safety first
set -eu
#set -eux

#echo "Script $0: count reads using grep '^@' | wc -l"  # OLD WAY
echo "Script $0: count reads reads by dividing number of lines by 4" 

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
	
	# Count reads in the zipped .fastq
	# OLD WAY: deprecated because fastq quality score line can start with @
	# append the number of @ lines in the file:
	# count_command="zcat $fq | grep '^@' | wc -l >> $RESULTS "
	# NEW WAY:
	count_command="echo $(zcat $fq | wc -l ) /4 | bc >> $RESULTS"
	echo "count command: $count_command"
	eval $count_command
	
done

echo "done counting all reads in each zipped fastq file"

source activate py3
set -eu

# Check for 2 arguments and help user if not
test $# -lt 2 && {
	echo "usage: $0 /work/rnaseq/rna_fastq_read_counts.tsv \"./*/*.sorted.bam\" bam_read_counts--may_over-report.tsv"
}

# /work/rnaseq/rna_fastq_read_counts.tsv
fastq_counts=$1

# get bam read counts
#bam_path='"./*/*.sorted.bam"'
bam_path=$2
#bam_counts=bam_read_counts--may_over-report.tsv
bam_counts=$3
bam_read_counter=/work/rnaseq/scripts/count_total_bam_reads.sh

# Look to see if the bams have already been counted.  Takes a while! (few hrs)
if [ -e "$bam_counts" ]
then
	echo "Found $bam_counts; don't re-compute"
else
	echo "Didn't find $bam_counts; compute now."
	echo 'Bam path to use: $bam_path'
	cmd="$bam_read_counter $bam_path $bam_counts"
	echo command: $cmd
	eval $cmd
fi

# Now use those bam counts to compare loss to fastq
python /work/rnaseq/scripts/analyze_read_loss--fastq_to_bam.py \
	$fastq_counts \
	$bam_counts \
	fastq_and_bam_reads.tsv



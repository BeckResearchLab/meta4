FASTQS='/work/rnaseq/data/*/Raw_Data/*.fastq.gz'
RESULTS=rna_fastq_read_counts.tsv
COUNT_SCRIPT=/work/m4b_binning/assembly/data/sample_info/count_reads_in_each_sample.sh

cmd="$COUNT_SCRIPT \"$FASTQS\" $RESULTS"
echo $cmd
eval $cmd



source deactivate
set -eu

read_mapping_script=/work/m4b_binning/assembly/map_reads_to_contigs/map_reads.sh
htseq_count=/work/m4b_binning/assembly/map_reads_to_contigs/htseq-count_reads.sh

reference_fasta=/work/m4b_binning/assembly/contigs_by_length/contigs_longer_than_1500bp.fa

#fastas_path='/work/rnaseq/data/Methane_oxidation_as_a_community_function__defining_partnerships_and_strategies_through_sequencing_metagenomes_and_metatranscriptomes_of_laboratory_manipulated_microcosms__*/Raw_Data/*.fastq.gz'
fastas_path='/work/rnaseq/data/simple_paths/*.fastq.gz'
merged_gff=/work/m4b_binning/assembly/prokka/contigs/contigs_longer_than_1500bp/contigs_longer_than_1500bp_gffs_concatenated.gff

fasta_files=`ls $fastas_path`

#echo "fasta files found:"
#for f in $fasta_files
#do
#	echo $f
#done

# Make the bam alignment file
#ls $fastas_path | parallel --dry-run --jobs 25 --results stdouterr_bwa_{} --joblog joblog_bwa.txt  $read_mapping_script {} $reference_fasta 
#ls $fastas_path | parallel --jobs 25 --results stdouterr_bwa_{} --joblog joblog_bwa.txt  $read_mapping_script {} $reference_fasta 
#ls $fastas_path | parallel --jobs 5 --results stdouterr_bwa_5jobs_{} --joblog joblog_bwa_5jobs.txt  $read_mapping_script {} $reference_fasta 
# 170217 do 4 jobs with 8 threads each, to reproduce what Dave did on the way to computing coverage.
# Also, map_reads.sh now asks for a threads argument.
ls $fastas_path | parallel --jobs 4 --results stdouterr_bwa_5jobs_{} --joblog joblog_bwa_5jobs.txt  $read_mapping_script {} $reference_fasta 8

bams_path=map_to_contigs_longer_than_1500bp/*/*.sorted.bam

#echo "bam files found:"
#for f in $bams_path
#do
#	echo $f
#done

# Run htseq-count
# ls $files | parallel --jobs 3 --results stdouterr_{} --joblog joblog.txt   ./toy_script.sh {} parallel_results.txt
#ls $bams_path | parallel --dry-run --jobs 25 --results stdouterr_htseq-count_{} --joblog joblog_htseq_count.txt  $htseq_count {} $merged_gff
#ls $bams_path | parallel --jobs 25 --results stdouterr_htseq-count_{} --joblog joblog_htseq_count.txt  $htseq_count {} $merged_gff
# Original (pre 2/17/2016):
#ls $bams_path | parallel --jobs 5 --results stdouterr_htseq-count_5jobs_{} --joblog joblog_htseq_count_5jobs.txt  $htseq_count {} $merged_gff
#Re-run by hand
ls $bams_path | parallel --jobs 5 --results stdouterr_htseq-count_5jobs_{} --joblog joblog_htseq_count_5jobs.txt  $htseq_count {} $merged_gff

# Aggregate 
python aggregate_counts.py map_to_contigs_longer_than_1500bp map_to_contigs_longer_than_1500bp.tsv

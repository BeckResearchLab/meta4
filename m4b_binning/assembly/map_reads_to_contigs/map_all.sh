set -eu

fastq_files=`ls /work/m4b_binning/assembly/data/*OW*/*.fastq.gz` # `` same as $()
echo "fastq files found:"
# You don't alwys need to echo something in quotes.
echo "$fastq_files" | wc -w

# Doesn't work: it isn't splititng up the fastq files. 
#echo $fastq_files | parallel --jobs 20 map_reads.sh {} /work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa

# test that worked:
#ls /work/m4b_binning/assembly/data/*OW*/*.fastq.gz | parallel --jobs 20 ./args.sh {} /work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa

# 170217 do 4 jobs with 8 threads each, to reproduce what Dave did on the way to computing coverage. 
ls /work/m4b_binning/assembly/data/*OW*/*.fastq.gz | parallel --jobs 4 --results stdout_stderr_170217_{} --joblog joblog_170217  ./map_reads.sh {} /work/m4b_binning/assembly/contigs_by_length/contigs_longer_than_1500bp.fa  8

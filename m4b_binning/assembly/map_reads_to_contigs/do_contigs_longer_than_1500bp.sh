# TODO: reference run_all.sh to get the bam files
#./map_all.sh  # ??

# TODO: make sure the aggregated gff is available

./htseq-count_reads--all_samples.sh map_to_contigs_longer_than_1500bp/ /work/m4b_binning/assembly/prokka/contigs/contigs_longer_than_1500bp/contigs_longer_than_1500bp_gffs_concatenated.gff

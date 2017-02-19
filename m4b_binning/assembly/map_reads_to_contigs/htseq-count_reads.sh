source activate py2

echo "use $@: use htseq-count on a gff to get a table of gene counts for one sample"

# safety first
set -eu 

# commad-line args:
test $# -lt 2 && {
   echo "usage: $0 bam_path gff_path"
   #echo "e.g. $0 /work/m4b_binning/assembly/data/HOW10/8777.1.112183.AATAGG.fastq.gz /work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa"
   exit
}

bam_path=$1
gff_path=$2
gff_file=$(basename $gff_path)
bam_file=$(basename $bam_path)
bam_base=$(echo $bam_file | sed 's#\.fa.*##') # cut off everything after .fa 
gff_base=$(echo $gff_file | sed 's#\.gff.*##') # cut off everything after .fa 
#echo "bam base: $bam_base"
#echo "gff base: $gff_base"


# htseq-count -m intersection-nonempty -s no -t CDS -i ID  -f bam map_to_contigs_longer_than_1500bp_test_scale/8777.4.112218.TGCCAT/8777.4.112218.TGCCAT.fastq.gz.sorted.bam /work/m4b_binning/assembly/prokka/contigs/contigs_longer_than_1500bp/contigs_longer_than_1500bp_gffs_concatenated.gff > ./map_to_contigs_longer_than_1500bp_test_scale/8777.4.112218.TGCCAT/test_CDS.summary.dat

#echo "bam file: $bam_file; bam path: $bam_path"
#echo "gff file: $gff_path"

outdir=$(dirname $bam_path)
#echo "outdir: $outdir"
outfile=$outdir/${bam_base}_to_${gff_base}.summary.dat
echo "save output to $outfile"

# First check that the expected file doesn't exist. 
if [[ -f $outfile ]] ; then
    echo "File $outfile already exists.  Don't re-run htseq-count."
    exit
fi

# Call htseq-count
cmd="htseq-count -m intersection-nonempty -s no -t CDS -i ID -f bam $bam_path $gff_path > $outfile"
echo "--- command to run ---"
echo $cmd
echo "--- ---"
eval $cmd


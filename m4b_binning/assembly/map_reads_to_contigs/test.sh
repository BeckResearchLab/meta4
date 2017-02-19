set -eu 

threads=$1
#fastq_path=/work/m4b_binning/assembly/data/HOW10/8777.1.112183.AAGCGA.fastq.gz
fastq_path=/work/m4b_binning/assembly/data/HOW10/8777.1.112183.AATAGG.fastq.gz
contigs_path=/work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa

# parse 8777.1.112183.AAGCGA type string from path like ...data/HOW10/8777.1.112183.AAGCGA.fastq.gz
sample=$(echo $fastq_path | sed 's#.*/\(.*\).fastq.gz#\1#') 
echo "parsed $sample from $fastq_path"

# parse the contig_details and put these mappings in there
contigs_base=$(echo $contigs_path| sed 's#.*/\(.*\).fa.*#\1#')
echo "parsed $contigs_base from $contigs_path"

dir="map_to_${contigs_base}_thread_tests"
mkdir -p $dir

dir=$dir/${sample}_${threads}_threads
mkdir -p $dir
cd $dir

t_file='timestamps.txt'
touch $t_file
echo "$(date)" > $t_file 

cp $fastq_path .
cp $contigs_path .
fastq=$(basename "$fastq_path")
contigs=$(basename "$contigs_path")
echo "do stuff on $fastq"

if [ ! -e $sample.bwt ]
then 
	bwa index $contigs
	echo "done with bwa index"
	samtools faidx $contigs
	echo "done with bwa samtools faidx"
fi

## the first part maps the fastqgz to the bin fasta that pipes to samtools to create a binary that is fed to another samtools to sort by contig
cmd="bwa mem -t $threads $contigs $fastq | samtools view -h -b -S /dev/stdin | samtools sort -m 1000000000 -o $fastq.sorted.bam /dev/stdin"
echo $cmd
eval $cmd

echo "$(date)" >> $t_file 

cd ..

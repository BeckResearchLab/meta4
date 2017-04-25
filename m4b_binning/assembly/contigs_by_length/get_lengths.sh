source activate py3
set -eu

SCRIPT=/work/general_scripts/make_length_tsv.py
FASTA=contigs_longer_than_1500bp.fa
python $SCRIPT $FASTA > $FASTA.lengths.tsv

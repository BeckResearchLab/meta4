source activate py3
set -eu

mkdir -p reuslts
mkdir -p results/data

script=run_GeneNet_and_summarize.R

input_fname=input_to_R--no_features_dropped.tsv
echo filename for all features as percent of fastq: $input_fname
python prep_input.py --dest_filename $input_fname > $script.log 2>&1

percent=0.001
Rscript $script -f $input_fname -p $percent >> $script.log 2>&1
# 170410 took half of the machine memory.  117681/245859MB --> 117GB
# This 30k input gene run seems noticeably slower than for the 28k version prior to the trimming fix



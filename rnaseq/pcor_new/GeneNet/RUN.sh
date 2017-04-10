source activate py3
set -eu
python prep_input.py

script=run_GeneNet_and_summarize.R
input=input_for_R--min_percent_0.005--unnormalized.tsv
percent=0.005
Rscript $script -f $input -p $percent > $script.log 2>&1

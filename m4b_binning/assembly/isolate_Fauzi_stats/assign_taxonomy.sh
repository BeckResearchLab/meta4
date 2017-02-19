# phylophlan needs BioPython
source activate py2

# need to refresh bashrc stuff
source ~/.bashrc

# safety first
set -eu

PHYLOPHLAN_DIR=/work/software/phylophlan/
#cd /work/software/phylophlan/nsegata-phylophlan-8e2d2ec74872
cd $PHYLOPHLAN_DIR

# phylophlan wants a folder in /phylophlan/input/my_project_name with .faa files.
# we have .fa files
mkdir -p ./input/fauz

# remove any old files:
rm -f ./Fauzi/bins/*.faa  # need `-f` when `set -eu`

for file in $(echo /work/m4b_binning/assembly/prokka/results/bin*/*.faa) 
do
 echo $file
 faa=${file##*/}  #  remove all dirs (same as "basename $FILE")
 echo "faa: " $faa
 dest_path=$PHYLOPHLAN_DIR/input/metabat/
 echo "dest path: "$dest_path
 faa_path=$dest_path$faa
 echo "faa_path: " $faa_path
 echo "cp $file $faa_path"

 # the names look good, so copy them over. 
 cp $file $faa_path
done 

#./phylophlan.py example_corynebacteria
echo "pwd: "$(pwd)
./phylophlan.py -i -t metabat

# move desired results into local folder 
mkdir -p /work/m4b_binning/assembly/phylophlan/
cp $PHYLOPHLAN_DIR/output/metabat/* /work/m4b_binning/assembly/phylophlan/

# ask it to make taxonomy
# ./phylophlan -i metabat
# *** WARNING *** Assuming DNA (see -seqtype option), invalid letters found: EFILPQ
#
#./phylophlan.py -t metabat

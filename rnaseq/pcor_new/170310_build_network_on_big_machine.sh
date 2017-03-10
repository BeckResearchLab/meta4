source activate py3

set -eu

output_dir=Ledoit_Wolf_large_instance_result
mkdir -p $output_dir
log=$output_dir/$output_dir.log
echo log output info to $log

#python network_prep.py -d $output_dir --list 1 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 5e-4 1e-4 5e-5 1e-5 1e-10  2>&1 | tee > $log
python network_prep.py -d $output_dir --list 1 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 5e-4 1e-4 5e-5 1e-5 1e-10  2>&1  > $log

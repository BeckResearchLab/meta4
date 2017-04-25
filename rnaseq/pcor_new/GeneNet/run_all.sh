set -eu

cd cutoff_0.1
./RUN.sh
cd ..
echo -----------------------------------

cd cutoff_0.05
./RUN.sh
cd ..
echo -----------------------------------

cd cutoff_0.01
./RUN.sh
cd ..
echo -----------------------------------

cd cutoff_0.005
./RUN.sh
cd ..
echo -----------------------------------

cd cutoff_0.001
./RUN.sh
cd ..
echo -----------------------------------

cd cutoff_0.0005
./RUN.sh
cd ..

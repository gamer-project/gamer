cd src
make
cd ../bin/Models
cp ../gamer gamer
sh clean.sh
export LD_LIBRARY_PATH=/home/return/miniconda3/pkgs/hdf5-1.10.4-hb1b8bf9_0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/return/gsl/lib:$LD_LIBRARY_PATH
mpirun --use-hwthread-cpus -np 1 -map-by ppr:1:socket:pe=4 ./gamer
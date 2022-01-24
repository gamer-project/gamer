if [ "$1" = "master" ];then
   git ck gamer-project/master
   MAKEFILE=MyMakefile_master
fi

if [ "$1" = "adaptive-minmod" ];then
   git ck adaptive-minmod
   MAKEFILE=MyMakefile
fi

# Compile GPU
if [ "$2" = "compile-y" ];then
   cd /projectW/tseng/gamer/adaptive-minmod-branch/src
   sed -i -e 's/#SIMU_OPTION += -DGPU/SIMU_OPTION += -DGPU/g' $MAKEFILE
   make clean
   rm /projectW/tseng/gamer/adaptive-minmod-branch/bin-gpu/gamer
   make -f $MAKEFILE -j 10
   cp /projectW/tseng/gamer/adaptive-minmod-branch/src/gamer /projectW/tseng/gamer/adaptive-minmod-branch/bin-gpu
fi

# Run GPU
cd /projectW/tseng/gamer/adaptive-minmod-branch/bin/fail-gpu
./clean.sh
mpirun  -map-by ppr:4:socket:pe=4 --bind-to core --report-bindings ./gamer >& log


# Compile CPU
if [ "$2" = "compile-y" ];then
   cd /projectW/tseng/gamer/adaptive-minmod-branch/src
   sed -i -e 's/SIMU_OPTION += -DGPU/#SIMU_OPTION += -DGPU/g' $MAKEFILE
   make clean
   rm /projectW/tseng/gamer/adaptive-minmod-branch/bin-cpu/gamer
   make -f $MAKEFILE -j 10
   cp /projectW/tseng/gamer/adaptive-minmod-branch/src/gamer /projectW/tseng/gamer/adaptive-minmod-branch/bin-cpu
fi

# Run CPU
cd /projectW/tseng/gamer/adaptive-minmod-branch/bin/fail-cpu
./clean.sh
mpirun  -map-by ppr:4:socket:pe=4 --bind-to core --report-bindings ./gamer >& log

GPU_FILE=/projectW/tseng/gamer/adaptive-minmod-branch/bin/fail-gpu/Data_000001
CPU_FILE=/projectW/tseng/gamer/adaptive-minmod-branch/bin/fail-cpu/Data_000001

GPU_OVER=`tail /projectW/tseng/gamer/adaptive-minmod-branch/bin/fail-gpu/Record__Dump -n 1 | awk '{print $1}'`
CPU_OVER=`tail /projectW/tseng/gamer/adaptive-minmod-branch/bin/fail-cpu/Record__Dump -n 1 | awk '{print $1}'`

while [[ $GPU_OVER != 1 ]] || [[ $CPU_OVER != 1 ]]
do
  GPU_OVER=`tail /projectW/tseng/gamer/adaptive-minmod-branch/bin/fail-gpu/Record__Dump -n 1 | awk '{print $1}'`
  CPU_OVER=`tail /projectW/tseng/gamer/adaptive-minmod-branch/bin/fail-cpu/Record__Dump -n 1 | awk '{print $1}'`
  sleep 1
  echo "Wait!!"
  date
done

echo "Done!!"

cd /projectW/tseng/gamer/adaptive-minmod-branch/tool/analysis/gamer_compare_data/Run
rm -rf log*
go.sh $CPU_FILE $GPU_FILE 1.0
go.sh $CPU_FILE $GPU_FILE 1e-3
go.sh $CPU_FILE $GPU_FILE 1e-7
go.sh $CPU_FILE $GPU_FILE 1e-10
go.sh $CPU_FILE $GPU_FILE 1e-16


# Plot
cd /projectW/tseng/gamer/adaptive-minmod-branch/bin/fail-cpu
python plot_scripts/Read__Parameters/general__ReadSliceMultiple.py -f plot__slice

echo "Done!!!!!!!!"

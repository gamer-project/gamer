set -e
NUM_THREADS=32 # number of threads of openMP
MY_HOST=`hostname`
MY_CLUSTER=${MY_HOST::(-2)}



#====================================================================================================
# Paths
#====================================================================================================
# eureka
if [[ ${MY_CLUSTER} != "eureka" ]]; then echo "ERROR : These paths are for eureka only!!"; exit 1; fi
export PATH=/software/hdf5/default/bin:$PATH
export PATH=/software/gsl/default/bin:$PATH
export PATH=/software/gcc/9.3.0/bin:$PATH
export PATH=/software/openmpi/4.1.1-ucx_mt-gnu-9.3.0/bin:$PATH

export LD_LIBRARY_PATH=/software/hdf5/default/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/projectW/tseng/opt/cfitsio/cfitsio-4.0.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/projectW/tseng/opt/CCfits/CCfits-2.6/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/projectW/tseng/opt/HEALPix/HEALPix-3.80/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/software/intel/oneapi/compiler/2021.1.1/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/projectW/tseng/opt/openmpi-gcc-9.3.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/software/gsl/default/lib:$LD_LIBRARY_PATH

# TW3
# module load compiler/intel/2020u4 OpenMPI/4.0.5



#====================================================================================================
# Compile
#====================================================================================================
make clean && make CFLAGS+=-DNUM_THREADS=${NUM_THREADS} CFLAGS+=-DXRAY_EROSITA
mv bin/Project bin/Project_XRay

make clean && make CFLAGS+=-DNUM_THREADS=${NUM_THREADS} CFLAGS+=-DSYNCHROTRON
mv bin/Project bin/Project_Synchrotron

make clean && make CFLAGS+=-DNUM_THREADS=${NUM_THREADS} CFLAGS+=-DLEPTONIC_GAMMARAY
mv bin/Project bin/Project_GammaRay

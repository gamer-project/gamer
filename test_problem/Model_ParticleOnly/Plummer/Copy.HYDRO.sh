#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_* End_*                        ../../../src/Init/
cp Flu_*                               ../../../src/Fluid/
cp Par_*                               ../../../src/Particle/
cp plot*.gpt                           ../../../bin/run/
cp Input__Flag_NParPatch               ../../../bin/run/
cp Input__TestProb                     ../../../bin/run/

cp Makefile.HYDRO                      ../../../src/Makefile
cp Input__Parameter.HYDRO              ../../../bin/run/Input__Parameter

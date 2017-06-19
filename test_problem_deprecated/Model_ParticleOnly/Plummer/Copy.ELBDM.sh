#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_TestProb.cpp End_TestProb.cpp  ../../../src/Init/
cp Flu_*                               ../../../src/Fluid/
cp Par_*                               ../../../src/Particle/
cp plot*.gpt                           ../../../bin/run/
cp Input__Flag_NParPatch               ../../../bin/run/
cp Input__TestProb                     ../../../bin/run/

cp Makefile.ELBDM                      ../../../src/Makefile
cp Input__Parameter.ELBDM              ../../../bin/run/Input__Parameter

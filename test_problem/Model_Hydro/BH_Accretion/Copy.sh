#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_TestProb.cpp                   ../../../src/Init/
cp Init_ExternalAcc.cpp                ../../../src/SelfGravity/
cp CUPOT_ExternalAcc.cu                ../../../src/SelfGravity/GPU_Gravity/
cp Flu_*                               ../../../src/Fluid/
cp Flag_UserCriteria.cpp               ../../../src/Refine/
cp End_*                               ../../../src/Init/
cp Aux_RecordUser.cpp                  ../../../src/Auxiliary/
cp CUFLU.h                             ../../../include/
cp Makefile                            ../../../src/
cp Input__*                            ../../../bin/run/

# make sure that CUPOT_ExternalAcc.1.1.0.cu will be recompiled
touch ../../../src/Model_Hydro/GPU_HydroGravity/CUPOT_HydroGravitySolver.cu

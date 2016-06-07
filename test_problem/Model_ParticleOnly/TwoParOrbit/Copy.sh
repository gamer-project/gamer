#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Makefile                            ../../../src/
cp Init_TestProb.cpp End_TestProb.cpp  ../../../src/Init/
cp Flu_*                               ../../../src/Fluid/
cp Par_*                               ../../../src/Particle/
cp Input__*                            ../../../bin/Run/
cp plot_*                              ../../../bin/Run/
cp Init_External*                      ../../../src/SelfGravity/
cp CUPOT_ExternalAcc.1.1.0.cu          ../../../src/SelfGravity/GPU_Gravity/
cp CUPOT_ExternalPot.1.0.0.cu          ../../../src/SelfGravity/GPU_Poisson/

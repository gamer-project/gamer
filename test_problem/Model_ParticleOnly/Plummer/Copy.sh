#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Makefile                            ../../../src/
cp Init_* End_*                        ../../../src/Init/
cp Flu_*                               ../../../src/Fluid/
cp Flag_*                              ../../../src/Refine/
cp Par_*                               ../../../src/Particle/
cp Input__*                            ../../../bin/Run/
cp plot*.gpt                           ../../../bin/Run/

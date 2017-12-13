#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_*                              ../../../src/Init/
cp Makefile                            ../../../src/
cp Input__*                            ../../../bin/Run/
cp Flu_Close.cpp                       ../../../src/Fluid/
cp ELBDM_Init_StartOver_AssignData.cpp ../../../src/Model_ELBDM/
cp Patch.h                             ../../../include/

#!/bin/bash

# copy files to the correct directories for the target test problem
cp Init_* End_*            ../../../src/Init/
cp Par_*                   ../../../src/Particle/
cp Flag_UserCriteria.cpp   ../../../src/Refine/
cp Makefile                ../../../src/Makefile
cp Input__*                ../../../bin/run/
cp *.dat                   ../../../bin/run/
cp plot*.py                ../../../bin/run/

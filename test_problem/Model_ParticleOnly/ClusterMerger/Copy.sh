#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_* End_*                        ../../../src/Init/
cp Par_*                               ../../../src/Particle/
cp plot*.gpt                           ../../../bin/run/
cp Input__*                            ../../../bin/run/
cp Makefile                            ../../../src/Makefile

#cp CUFLU.h                             ../../../include/

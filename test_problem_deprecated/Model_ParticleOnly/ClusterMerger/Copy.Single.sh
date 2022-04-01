#!/bin/bash

# copy files to the correct directories for the target test problem
cp Init_* End_*                        ../../../src/Init/
cp Par_*                               ../../../src/Particle/
cp plot*.gpt                           ../../../bin/run/
cp Makefile                            ../../../src/Makefile
cp Input__Flag_NParPatch               ../../../bin/run/

cp Input__Parameter.Single             ../../../bin/run/Input__Parameter
cp Input__TestProb.Single              ../../../bin/run/Input__TestProb

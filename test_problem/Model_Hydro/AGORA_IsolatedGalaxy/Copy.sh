#!/bin/bash

# copy files to the correct directories for the target test problem
cp Init_* End_*   ../../../src/Init/
cp Par_*          ../../../src/Particle/
cp Makefile       ../../../src/Makefile
cp Input__*       ../../../bin/run/
cp vcirc.dat      ../../../bin/run/

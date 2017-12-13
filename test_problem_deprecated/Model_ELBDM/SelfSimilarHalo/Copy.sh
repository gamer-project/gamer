#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_* End_*                        ../../../src/Init/
cp Makefile                            ../../../src/Makefile
cp Input__*                            ../../../bin/run/
cp plot_slice.py                       ../../../bin/run/

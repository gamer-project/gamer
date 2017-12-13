#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_*                              ../../../src/Init/
cp End_*                               ../../../src/Init/
cp Makefile                            ../../../src/
cp Flu_*                               ../../../src/Fluid/
cp Input__*                            ../../../bin/Run/ -L
cp plot.gpt                            ../../../bin/Run/

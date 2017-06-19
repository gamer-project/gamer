#!/bin/bash

# copy files to the correct directories for the target test problem
cp Init_TestProb.cpp                   ../../../src/Init/
cp End_*                               ../../../src/Init/
cp Makefile                            ../../../src/
cp Input__*                            ../../../bin/run/
cp plot*                               ../../../bin/run/

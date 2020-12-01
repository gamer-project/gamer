#!/bin/bash
make clean
make -f MyMakefile -j $1

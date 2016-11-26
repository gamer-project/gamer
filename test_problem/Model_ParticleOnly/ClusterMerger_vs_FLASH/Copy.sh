#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_* End_*   ../../../src/Init/
cp Par_*          ../../../src/Particle/
cp plot*.gpt      ../../../bin/run/
cp Makefile       ../../../src/Makefile
cp Input__*       ../../../bin/run/
cp CUFLU.h        ../../../include/

ln -s /home/hyschive/project/amr_code_comparison/merging_cluster/gamer/tool/convert_particle_hdf5-to-text/particle.cbin.f8 \
      ../../../bin/run/
ln -s /home/hyschive/project/amr_code_comparison/merging_cluster/gamer/tool/convert_profile_hdf5-to-text/profile.txt \
      ../../../bin/run/

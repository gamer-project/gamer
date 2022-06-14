#!/bin/bash

ffmpeg -f image2 -framerate 10 -pattern_type glob -i "Data_*_Slice_z_velocity_magnitude.png" -y -vcodec libx264 -crf 10  -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" velm.mp4

ffmpeg -f image2 -framerate 10 -pattern_type glob -i "Data_*_Slice_z_particle_density_on_grid.png" -y -vcodec libx264 -crf 10  -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" pden.mp4

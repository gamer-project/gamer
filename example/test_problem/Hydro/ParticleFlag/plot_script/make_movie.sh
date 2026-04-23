#!/bin/bash

ffmpeg -f image2 -framerate 2 -pattern_type glob -i "Data_*_Slice_z_density.png" -y -vcodec libx264 -crf 10 \
       -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" movie_par_flag.mp4

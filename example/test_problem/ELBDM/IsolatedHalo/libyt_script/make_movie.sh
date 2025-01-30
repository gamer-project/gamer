#!/usr/bin/bash

ffmpeg -framerate 5/1 -i Fig%06d_Slice_x_Dens_mode_c_central_halo.png -c:v libx264 -preset slow -tune animation -pix_fmt yuv420p -s 1920x1200 Slice_x_Dens_mode_central_halo.mp4

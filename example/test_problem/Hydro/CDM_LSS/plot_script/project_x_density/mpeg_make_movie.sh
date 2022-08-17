ffmpeg -framerate 5/1 -i Data_%06d_Particle_x_particle_mass.png -c:v libx264 -preset slow -tune animation \
       -pix_fmt yuv420p -s 1920x1200 CDM_projected_x_particle_density.mp4
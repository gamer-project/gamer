for N in gas_dens gas_temp; do

   mkdir $N
   cd $N
   mkdir slice-x slice-z proj-x proj-z
   cd ..

done

mkdir gas_vabs
mkdir gas_vabs/slice-x gas_vabs/slice-z
mkdir par_mass
mkdir par_mass/proj-x par_mass/proj-z
mkdir gas_profile
mkdir gas_temp-vs-dens


mv *Slice_x_density.png                   gas_dens/slice-x
mv *Slice_z_density.png                   gas_dens/slice-z
mv *Projection_x_density.png              gas_dens/proj-x
mv *Projection_z_density.png              gas_dens/proj-z

mv *Slice_x_temperature*.png              gas_temp/slice-x
mv *Slice_z_temperature*.png              gas_temp/slice-z
mv *Projection_x_temperature*.png         gas_temp/proj-x
mv *Projection_z_temperature*.png         gas_temp/proj-z

mv *Slice_x_velocity_magnitude.png        gas_vabs/slice-x
mv *Slice_z_velocity_magnitude.png        gas_vabs/slice-z

mv *Particle_x_*.png                      par_mass/proj-x
mv *Particle_z_*.png                      par_mass/proj-z

mv fig__gas_profile* gas_profile/

mv *density_temperature_cell_mass* gas_temp-vs-dens


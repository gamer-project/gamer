for N in gas_dens gas_temp; do

   mkdir $N
   cd $N
   mkdir slice-x slice-z proj-x proj-z
   cd ..

done

mkdir gas_vabs
mkdir gas_vabs/slice-x gas_vabs/slice-z
mkdir par_mass_all
mkdir par_mass_all/proj-x par_mass_all/proj-z
mkdir par_mass_new
mkdir par_mass_new/proj-x par_mass_new/proj-z
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

mv *AllPar*Particle_x_*.png               par_mass_all/proj-x
mv *AllPar*Particle_z_*.png               par_mass_all/proj-z

mv *NewPar*Particle_x_*.png               par_mass_new/proj-x
mv *NewPar*Particle_z_*.png               par_mass_new/proj-z

mv fig__gas_profile* gas_profile/

mv *density_temperature_cell_mass* gas_temp-vs-dens

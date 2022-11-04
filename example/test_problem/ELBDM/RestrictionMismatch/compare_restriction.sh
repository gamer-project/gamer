#Run three times to generate three runs (without refinement, with restriction, with new restriction)
#Plot differences between the three runs
cp Input__Parameter_NoRefinement Input__Parameter
./gamer
cp Data_000000 RESTART
cp Data_000001 Data_000001_NoRefinement
cp Data_000002 Data_000002_NoRefinement
cp Input__Parameter_PhaseResOff Input__Parameter
./gamer
cp Data_000001 Data_000001_OldRestriction
cp Data_000002 Data_000002_OldRestriction
cp Input__Parameter_PhaseResOn Input__Parameter
./gamer
cp Data_000001 Data_000001_NewRestriction
cp Data_000002 Data_000002_NewRestriction
python plot_comparison.py
python plot_slice.py -s 0 -e 2
rm Data_000000_Slice_z_density_gradient_magnitude.png
rm Data_000000_Slice_z_Phase.png
rm Data_000001_Slice_z_density_gradient_magnitude.png
rm Data_000001_Slice_z_Phase.png
rm Data_000002_Slice_z_density_gradient_magnitude.png
rm Data_000002_Slice_z_Phase.png

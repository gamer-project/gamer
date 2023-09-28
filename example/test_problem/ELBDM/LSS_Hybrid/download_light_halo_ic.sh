filename=Music_InitCondition_z99_L2.8_N0064_LightHalo
link=http://use.yt/upload/a67d8dd1

rm UM_IC*
curl ${link} -o ${filename}
python3 elbdm_rescale_periodic_IC.py -input ${filename} -output Music_InitCondition_z99_L2.8_N0256_LightHalo -n_in 64 -n_out 256
ln -sf Music_InitCondition_z99_L2.8_N0256_LightHalo UM_IC_wave_light
python3 elbdm_wave_to_hybrid_IC.py -input UM_IC_wave_light -output UM_IC_hybrid_light -resolution 256
ln -sf UM_IC_hybrid_light UM_IC

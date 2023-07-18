filename=Music_InitCondition_z99_L2.8_N0064_LightHalo
link=http://use.yt/upload/a67d8dd1

rm UM_IC*
curl ${link} -o ${filename}
python3 ELBDM_RESCALE_PERIODIC_IC.py -input ${filename} -output Music_InitCondition_z99_L2.8_N0256_LightHalo -n_in 64 -n_out 256
ln -s Music_InitCondition_z99_L2.8_N0256_LightHalo UM_IC_wave
python3 ELBDM_WAVE2HYBRID_IC.py -input UM_IC_wave -output UM_IC_hybrid -resolution 256
ln -s UM_IC_hybrid UM_IC

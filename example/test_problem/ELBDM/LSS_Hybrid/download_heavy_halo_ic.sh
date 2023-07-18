filename=Music_InitCondition_z99_L2.8_N0256_HeavyHalo
link=http://use.yt/upload/4587e2e6

curl ${link} -o ${filename}
ln -s ${filename} UM_IC_wave
python ELBDM_WAVE2HYBRID_IC.py -input UM_IC_wave -output UM_IC_hybrid -resolution 256
ln -s UM_IC_hybrid UM_IC

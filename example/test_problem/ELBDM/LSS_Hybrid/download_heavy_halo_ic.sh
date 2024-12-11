filename=Music_InitCondition_z99_L2.8_N0256_HeavyHalo
link=https://use.yt/upload/4587e2e6

curl -L ${link} -o ${filename}
ln -sf ${filename} UM_IC_wave_heavy
python3 elbdm_wave_to_hybrid_IC.py -input UM_IC_wave_heavy -output UM_IC_hybrid_heavy -resolution 256
ln -sf UM_IC_hybrid_heavy UM_IC

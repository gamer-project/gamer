filename=UM_IC_z100_L2.8_N256_seed_1578581
link=https://girder.hub.yt/api/v1/item/675fee70999605c485c8de12/download

curl ${link} -o ${filename}
ln -sf ${filename} UM_IC_wave
python3 elbdm_wave_to_hybrid_IC.py -input UM_IC_wave -output UM_IC_hybrid -resolution 256
ln -sf UM_IC_hybrid UM_IC

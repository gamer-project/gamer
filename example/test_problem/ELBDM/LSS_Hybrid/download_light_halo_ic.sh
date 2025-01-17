#!/bin/bash

LOCAL_FILENAME="Music_InitCondition_z99_L2.8_N0064_LightHalo"
FILE_ID="6780d97e999605c485c8dea6"

# 1. clean
rm UM_IC*

# 2. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 3. link
python3 elbdm_rescale_periodic_IC.py -input ${LOCAL_FILENAME} -output Music_InitCondition_z99_L2.8_N0256_LightHalo -n_in 64 -n_out 256
ln -sf Music_InitCondition_z99_L2.8_N0256_LightHalo UM_IC_wave_light
python3 elbdm_wave_to_hybrid_IC.py -input UM_IC_wave_light -output UM_IC_hybrid_light -resolution 256
ln -sf UM_IC_hybrid_light UM_IC

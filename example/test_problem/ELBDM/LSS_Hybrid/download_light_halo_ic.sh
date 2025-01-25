#!/bin/bash

LOCAL_FILENAME="Music_InitCondition_z100_L2.8_N0064_LightHalo"
FILE_ID="6780d97e999605c485c8dea6"
FILE_SHA256="89e1626405e38e3bd756167e3c06719f193456aa8e7d210ddfd13246df0e1a91"

# 1. clean
rm UM_IC*

# 2. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 3. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"

# 4. link
python3 elbdm_rescale_periodic_IC.py -input ${LOCAL_FILENAME} -output Music_InitCondition_z100_L2.8_N0256_LightHalo -n_in 64 -n_out 256
ln -sf Music_InitCondition_z100_L2.8_N0256_LightHalo UM_IC_wave_light
python3 elbdm_wave_to_hybrid_IC.py -input UM_IC_wave_light -output UM_IC_hybrid_light -resolution 256
ln -sf UM_IC_hybrid_light UM_IC

#!/bin/bash

LOCAL_FILENAME="Music_InitCondition_z100_L2.8_N0256_HeavyHalo"
FILE_ID="6780d9d2999605c485c8dea9"
FILE_SHA256="2f15920763e6189abd81b6f39fd283ebeeb6b6b90dbbce0eb898d098ef4b497d"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"

# 3. link
ln -sf ${LOCAL_FILENAME} UM_IC_wave_heavy
python3 elbdm_wave_to_hybrid_IC.py -input UM_IC_wave_heavy -output UM_IC_hybrid_heavy -resolution 256
ln -sf UM_IC_hybrid_heavy UM_IC

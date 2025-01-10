#!/bin/bash

LOCAL_FILENAME="Music_InitCondition_z3200_L1.4_N0256_s1002"
FILE_ID="6780d8d6999605c485c8dea0"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 2. link
ln -s ${LOCAL_FILENAME} UM_IC

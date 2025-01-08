#!/bin/bash

LOCAL_FILENAME="HALO_IC_m22_1_Mh_4e9"
FILE_ID="677cc8db999605c485c8de83"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 2. link
ln -s ${LOCAL_FILENAME} HALO_IC_Halo1
ln -s ${LOCAL_FILENAME} HALO_IC_Halo2

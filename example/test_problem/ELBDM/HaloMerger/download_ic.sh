#!/bin/bash

LOCAL_FILENAME="HALO_IC_m22_1_Mh_4e9"
FILE_ID="677cc8db999605c485c8de83"
FILE_SHA256="acbd85842de65ff2360c7f3a1d1101c6f4f8939f430c3f61b8bc5f6f9a72fe94"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"

# 3. link
ln -s ${LOCAL_FILENAME} HALO_IC_Halo1
ln -s ${LOCAL_FILENAME} HALO_IC_Halo2

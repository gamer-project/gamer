#!/bin/bash

LOCAL_FILENAME="disk-heating-ic"
FILE_ID="677dd2d0999605c485c8de8f"
FILE_SHA256="5c981ffe1f0cd85237b51e9e2872e8047dad8a87e0419575255e4c1d5d8cf17a"

# 1 download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}.tgz"

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME}.tgz | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}.tgz"

# 3. unzip and link
tar -zxvf ${LOCAL_FILENAME}.tgz
rm ${LOCAL_FILENAME}.tgz
ln -s ${LOCAL_FILENAME}/UM_IC_0.4_M7 UM_IC
ln -s ${LOCAL_FILENAME}/PAR_IC_0.4_M7_low_res DiskHeatingParticleIC

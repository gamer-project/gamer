#!/bin/bash

LOCAL_FILENAME="disk-heating-ic"
FILE_ID="677dd2d0999605c485c8de8f"

# 1 download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}.tgz"

# 2. unzip and link
tar -zxvf ${LOCAL_FILENAME}.tgz
rm ${LOCAL_FILENAME}.tgz
ln -s ${LOCAL_FILENAME}/UM_IC_0.4_M7 UM_IC
ln -s ${LOCAL_FILENAME}/PAR_IC_0.4_M7_low_res PAR_IC

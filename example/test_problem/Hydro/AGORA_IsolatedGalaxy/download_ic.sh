#!/bin/bash

LOCAL_FILENAME1="LOW"
LOCAL_FILENAME2="CloudyData_UVB=HM2012.h5"
FILE_ID1="677ca225999605c485c8de6f"
FILE_ID2="677ca211999605c485c8de6c"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID1}/download -o "${LOCAL_FILENAME1}.tar.gz"
curl https://hub.yt/api/v1/item/${FILE_ID2}/download -o "${LOCAL_FILENAME2}"

# 2. unzip
tar xzvf ${LOCAL_FILENAME1}.tar.gz
mv ${LOCAL_FILENAME1}/*.dat ./
rmdir ${LOCAL_FILENAME1}
rm ${LOCAL_FILENAME1}.tar.gz

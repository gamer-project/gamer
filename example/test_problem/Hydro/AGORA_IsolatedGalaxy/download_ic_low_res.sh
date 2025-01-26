#!/bin/bash

LOCAL_FILENAME1="LOW"
LOCAL_FILENAME2="CloudyData_UVB=HM2012.h5"
FILE_ID1="677ca225999605c485c8de6f"
FILE_ID2="677ca211999605c485c8de6c"
FILE_SHA256_1="8fd18d016f5af665d31d853a394c9089503df225f4c256eb1163c9a7b4e94644"
FILE_SHA256_2="8715f1b39e90a7296ec2adcd442fa13a3d45d2ad021c6fa2fae9e4ab7a4700b2"

# file download
curl https://hub.yt/api/v1/item/${FILE_ID1}/download -o "${LOCAL_FILENAME1}.tar.gz"
curl https://hub.yt/api/v1/item/${FILE_ID2}/download -o "${LOCAL_FILENAME2}"

# compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME1}.tar.gz | awk '{print $1}'` = "${FILE_SHA256_1}" ] && echo "File broken: ${LOCAL_FILENAME1}.tar.gz"
! [ `sha256sum ${LOCAL_FILENAME2}        | awk '{print $1}'` = "${FILE_SHA256_2}" ] && echo "File broken: ${LOCAL_FILENAME2}"

# file unzip
tar xzvf ${LOCAL_FILENAME1}.tar.gz
mv ${LOCAL_FILENAME1}/*.dat ./
rmdir ${LOCAL_FILENAME1}
rm ${LOCAL_FILENAME1}.tar.gz

# Input_* soft links
ln -fs ./Input_Options/Input__Flag_Jeans.low-res Input__Flag_Jeans
ln -fs ./Input_Options/Input__Flag_ParMassCell.low-res Input__Flag_ParMassCell
ln -fs ./Input_Options/Input__Flag_Rho.low-res Input__Flag_Rho
ln -fs ./Input_Options/Input__Parameter.low-res Input__Parameter

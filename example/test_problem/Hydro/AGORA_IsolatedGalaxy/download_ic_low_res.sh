#!/bin/bash

LOCAL_FILENAME1="LOW"
LOCAL_FILENAME2="CloudyData_UVB=HM2012.h5"
FILE_ID1="677ca225999605c485c8de6f"
FILE_ID2="677ca211999605c485c8de6c"

# file download
curl https://hub.yt/api/v1/item/${FILE_ID1}/download -o "${LOCAL_FILENAME1}.tar.gz"
curl https://hub.yt/api/v1/item/${FILE_ID2}/download -o "${LOCAL_FILENAME2}"

# file unzip
tar xzvf ${LOCAL_FILENAME1}.tar.gz
mv ${LOCAL_FILENAME1}/*.dat ./
rmdir ${LOCAL_FILENAME1}
rm ${LOCAL_FILENAME1}.tar.gz

# Input_* soft links
ln -s ./Input_Options/Input__Flag_Jeans.low-res Input__Flag_Jeans
ln -s ./Input_Options/Input__Flag_ParMassCell.low-res Input__Flag_ParMassCell
ln -s ./Input_Options/Input__Flag_Rho.low-res Input__Flag_Rho
ln -s ./Input_Options/Input__Parameter.low-res Input__Parameter

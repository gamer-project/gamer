#!/bin/bash
LOCAL_FILENAME1="HI"
LOCAL_FILENAME2="CloudyData_UVB=HM2012.h5"
FILE_ID1="677e4757999605c485c8de92"
FILE_ID2="677ca211999605c485c8de6c"

# file download
curl https://hub.yt/api/v1/item/${FILE_ID1}/download -o "${LOCAL_FILENAME1}.tar.gz"
curl https://hub.yt/api/v1/item/${FILE_ID2}/download -o "${LOCAL_FILENAME2}"

# file unzip
tar xzvf ${LOCAL_FILENAME1}.tar.gz
mv ${LOCAL_FILENAME1}/* ./
rm -rf ${LOCAL_FILENAME1}
rm ${LOCAL_FILENAME1}.tar.gz

# Input_* soft links
ln -s ./Input_Options/Input__Flag_Jeans.high-res Input__Flag_Jeans
ln -s ./Input_Options/Input__Flag_ParMassCell.high-res Input__Flag_ParMassCell
ln -s ./Input_Options/Input__Flag_Rho.high-res Input__Flag_Rho
ln -s ./Input_Options/Input__Parameter.high-res Input__Parameter

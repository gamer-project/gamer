#!/bin/bash
# file download
wget --no-check-certificate -O ./LOW.tar.gz https://www.dropbox.com/sh/1xzt1rysy9v3a9l/AAAMlJBQG1OQFW4cjhp11Ex6a/LOW.tar.gz?dl=1
wget --no-check-certificate  -O ./CloudyData_UVB=HM2012.h5 https://github.com/grackle-project/grackle_data_files/raw/main/input/CloudyData_UVB=HM2012.h5

# file unzip
tar xzvf LOW.tar.gz
mv LOW/*.dat ./
rmdir LOW
rm LOW.tar.gz

# Input_* soft links
ln -s ./Input_Options/Input__Flag_Jeans.low-res Input__Flag_Jeans
ln -s ./Input_Options/Input__Flag_ParMassCell.low-res Input__Flag_ParMassCell
ln -s ./Input_Options/Input__Flag_Rho.low-res Input__Flag_Rho
ln -s ./Input_Options/Input__Parameter.low-res Input__Parameter

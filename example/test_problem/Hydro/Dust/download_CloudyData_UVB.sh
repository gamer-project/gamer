#!/bin/bash
LOCAL_FILENAME="CloudyData_UVB=HM2012.h5"

# file download
curl -L https://raw.githubusercontent.com/grackle-project/grackle_data_files/main/input/CloudyData_UVB=HM2013.h5 -o "${LOCAL_FILENAME}"

#!/bin/bash

API_URL="https://girder.hub.yt/api/v1"
FILE_ID="66ff5e05999605c485c8d67a"
LOCAL_FILE="Zoomin_IC"

# download
girder-cli --api-url ${API_URL} download --parent-type item ${FILE_ID} ${LOCAL_FILE}

# unzip
tar zxvf ${LOCAL_FILE}/Zoomin_IC_heavyhalo.tar.gz
rm -r ${LOCAL_FILE}
ln -fs UM_IC_hybrid_N1024 UM_IC

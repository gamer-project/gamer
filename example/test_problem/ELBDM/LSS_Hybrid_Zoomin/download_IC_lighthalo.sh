#!/bin/bash

API_URL="https://girder.hub.yt/api/v1"
FILE_ID="66ff5e15999605c485c8d67d"
LOCAL_FILE="Zoomin_IC"

# download
girder-cli --api-url ${API_URL} download --parent-type item ${FILE_ID} ${LOCAL_FILE}

# unzip
tar zxvf ${LOCAL_FILE}/Zoomin_IC_lighthalo.tar.gz
rm -r ${LOCAL_FILE}
ln -fs UM_IC_hybrid_N256_third UM_IC

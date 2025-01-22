#!/bin/bash

LOCAL_FILENAME="UM_IC_run05-halo08-lv4"
FILE_ID="677cbad6999605c485c8de77"
FILE_SHA256="7ed91ba48a9aec139e0574629b689090ae43496fb957c6822c7ec1bd1217e22e"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}.tgz"

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME}.tgz | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}.tgz"

# 3. unzip and link
tar -zxvf ${LOCAL_FILENAME}.tgz
rm ${LOCAL_FILENAME}.tgz
ln -s ${LOCAL_FILENAME} UM_IC

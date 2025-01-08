#!/bin/bash

LOCAL_FILENAME="UM_IC_run05-halo08-lv4"
FILE_ID="677cbad6999605c485c8de77"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}.tgz"

# 2. unzip and link
tar -zxvf ${LOCAL_FILENAME}.tgz
rm ${LOCAL_FILENAME}.tgz
ln -s ${LOCAL_FILENAME} UM_IC

#!/bin/bash

LOCAL_FILENAME="uniform-granule"
FILE_ID="67da7961ef64ad0f8e84e795"

# 1. download
curl -L https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}.tar.gz"

# 2. unzip
tar -zxvf ${LOCAL_FILENAME}.tar.gz
rm ${LOCAL_FILENAME}.tar.gz


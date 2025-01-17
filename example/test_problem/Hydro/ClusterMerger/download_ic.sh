#!/bin/bash

LOCAL_FILENAME="gamer_ic_merging_cluster.tgz"
FILE_ID="677caaec999605c485c8de74"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}.tar.gz"

# 2. unzip
tar -zxvf ${LOCAL_FILENAME}
rm ${LOCAL_FILENAME}

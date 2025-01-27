#!/bin/bash

LOCAL_FILENAME="gamer_ic_merging_cluster.tgz"
FILE_ID="677caaec999605c485c8de74"
FILE_SHA256="a233a892818504cf15e188bca862e22250bb1f3e09155740e45d272e4ab5f1c1"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"

# 3. unzip
tar -zxvf ${LOCAL_FILENAME}
rm ${LOCAL_FILENAME}

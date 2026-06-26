#!/bin/bash

LOCAL_FILENAME="gamer_ic_cool_core_destruction.tgz"
FILE_ID="693f97b2467de92389bb7d15"
FILE_SHA256="138d5a6f594f4f3051b1c75f5f2837b4fe272790794766d80770653d58c33d45"

# 1. download
curl -L https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"

# 3. unzip
tar -zxvf ${LOCAL_FILENAME}
rm ${LOCAL_FILENAME}

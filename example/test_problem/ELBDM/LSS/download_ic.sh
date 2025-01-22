#!/bin/bash

LOCAL_FILENAME="Music_InitCondition_z3200_L1.4_N0256_s1002"
FILE_ID="6780d8d6999605c485c8dea0"
FILE_SHA256="114fd2a0d37e70ba7bd06907c878bd1c752ea76882ad1833af165696eef8cf9d"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"

# 3. link
ln -s ${LOCAL_FILENAME} UM_IC

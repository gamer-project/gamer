#!/bin/bash
# produced by MUSIC with GAMER plug-in (MUSIC commit 6248d133ab20aa84c9cb3b843a8e021a3119be70)

LOCAL_FILENAME="PAR_IC"
FILE_ID="6777e0e8999605c485c8de4b"
FILE_SHA256="80549dc26c6896b0f1efb76a378aec05d5fb983b22d6fd0c5323862e5ed2b7b9"

# 1. download
curl -L https://hub.yt/api/v1/item/${FILE_ID}/download -o ${LOCAL_FILENAME}

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"
#!/bin/bash
# produced by MUSIC with GAMER plug-in (MUSIC commit 6248d133ab20aa84c9cb3b843a8e021a3119be70)

LOCAL_FILENAME="PAR_IC"
FILE_ID="69e9d7524a11f772ce54c597"
FILE_SHA256="1c06000f3ae1d889434ec4d2f28fcb502dd73493eb28e0ebca2e54d072128afb"

# 1. download
curl -L https://girder.hub.yt/api/v1/file/${FILE_ID}/download -o ${LOCAL_FILENAME}

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"
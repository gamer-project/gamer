#!/bin/bash

LOCAL_FILENAME="spectral_tables.zip"
FILE_ID="6780d950999605c485c8dea3"
FILE_SHA256="304fb4d098d6ad6f6533f137fc78a4d05d2abf7c239392be29f694503410247f"

# 1. clean
rm -r spectral_tables*

# 2. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 3. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"

# 4. unzip
unzip ${LOCAL_FILENAME}
rm ${LOCAL_FILENAME}

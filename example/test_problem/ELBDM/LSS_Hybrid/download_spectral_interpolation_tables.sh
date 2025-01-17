#!/bin/bash

LOCAL_FILENAME="spectral_tables.zip"
FILE_ID="6780d950999605c485c8dea3"

# 1. clean
rm -r spectral_tables*

# 2. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# 3. unzip
unzip ${LOCAL_FILENAME}

#!/bin/bash

LOCAL_FILENAME="PAR_IC"
FILE_ID="677c92db999605c485c8de69"

# 1. Download through `curl`
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o ${LOCAL_FILENAME}

# 2. download through `girder`
API_URL="https://girder.hub.yt/api/v1"

# girder-cli --api-url ${API_URL} download --parent-type item ${FILE_ID} temp
# mv temp/${LOCAL_FILENAME} ./
# rmdir temp

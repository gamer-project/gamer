#!/bin/bash

LOCAL_FILENAME="PAR_IC"

# 1. Download through `curl`
curl https://hub.yt/api/v1/item/677c92db999605c485c8de69/download -o ${LOCAL_FILENAME}

# 2. download through `girder`
API_URL="https://girder.hub.yt/api/v1"
FILE_ID="677ba3e4999605c485c8de5a"

# girder-cli --api-url ${API_URL} download --parent-type item ${FILE_ID} ${LOCAL_FILENAME}


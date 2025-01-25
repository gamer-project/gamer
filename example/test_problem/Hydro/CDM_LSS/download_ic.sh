#!/bin/bash

LOCAL_FILENAME="PAR_IC"
FILE_ID="677c92db999605c485c8de69"
FILE_SHA256="46e27324953bcd7b4eecaecd395b6cf6ccbf662e65a12c935b3823abd8119be3"

# 1. download
curl https://hub.yt/api/v1/item/${FILE_ID}/download -o ${LOCAL_FILENAME}

# 2. compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME} | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"

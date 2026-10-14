#!/bin/bash
LOCAL_FILENAME="CloudyData_UVB=HM2012.h5"
FILE_ID="677ca211999605c485c8de6c"
FILE_SHA256="8715f1b39e90a7296ec2adcd442fa13a3d45d2ad021c6fa2fae9e4ab7a4700b2"

# file download
curl -L https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}"

# compare sha256sum
! [ `sha256sum ${LOCAL_FILENAME}        | awk '{print $1}'` = "${FILE_SHA256}" ] && echo "File broken: ${LOCAL_FILENAME}"

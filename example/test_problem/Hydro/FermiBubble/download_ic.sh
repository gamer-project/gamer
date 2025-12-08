#!/bin/bash

LOCAL_FILENAME="FermiBubble"
FILE_ID="677d3e9e999605c485c8de8c"

# 1. download
curl -L https://hub.yt/api/v1/item/${FILE_ID}/download -o "${LOCAL_FILENAME}.tar.gz"

# 2. unzip
tar zxvf ${LOCAL_FILENAME}.tar.gz
rm ${LOCAL_FILENAME}.tar.gz
mv IC/FermiBubble_IC ./
mv IC/R12 ./
rm -rf IC

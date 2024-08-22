API_URL="https://girder.hub.yt/api/v1"
FILE_ID="66c6bdb532f323dee1b80149"
LOCAL_FOLDER="./"

# download
girder-cli --api-url ${API_URL} download --parent-type item ${FILE_ID} ${LOCAL_FOLDER}

# unzip
tar zxvf FermiBubble.tar.gz
rm FermiBubble.tar.gz
mv IC/FermiBubble_IC ./
mv IC/R12 ./
rm -rf IC

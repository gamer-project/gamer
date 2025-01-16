filename=disk-heating-ic

curl https://girder.hub.yt/api/v1/item/6645cffcff473673ea91b24d/download -o ${filename}.tgz
tar -zxvf ${filename}.tgz
rm ${filename}.tgz
ln -s ${filename}/UM_IC_0.4_M7 UM_IC
ln -s ${filename}/PAR_IC_0.4_M7_low_res DiskHeatingParticleIC



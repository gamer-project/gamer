filename=ic_files

curl http://use.yt/upload/edc16152 -o ${filename}.tgz
tar -zxvf ${filename}.tgz
rm ${filename}.tgz
ln -s ${filename}/UM_IC_0.4_M7 UM_IC

curl http://use.yt/upload/3dde6b5c -o ${filename}.tgz
tar -zxvf ${filename}.tgz -C ${filename}/
rm ${filename}.tgz
ln -s ${filename}/PAR_IC_0.4_M7_low_res PAR_IC



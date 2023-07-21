filename=ic_files

curl http://use.yt/upload/edc16152 -o ${filename}.tgz 
tar -zxvf ${filename}.tgz
rm ${filename}.tgz

ln -s ${filename}/UM_IC_0.4_M7 UM_IC 
ln -s ${filename}/PAR_IC_0.4_M7 PAR_IC 



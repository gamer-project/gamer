filename=UM_IC_run05-halo08-lv4
link=http://use.yt/upload/d75aa595

curl ${link} -o ${filename}.tgz
tar -zxvf ${filename}.tgz
rm ${filename}.tgz

ln -s ${filename} UM_IC

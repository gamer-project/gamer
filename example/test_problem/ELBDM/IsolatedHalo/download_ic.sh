filename=UM_IC_run05-halo08-lv4
link=https://use.yt/upload/d75aa595

curl -L ${link} -o ${filename}.tgz
tar -zxvf ${filename}.tgz
rm ${filename}.tgz

ln -s ${filename} UM_IC

wget -O ./LOW.tar.gz https://www.dropbox.com/sh/1xzt1rysy9v3a9l/AAAMlJBQG1OQFW4cjhp11Ex6a/LOW.tar.gz?dl=1
wget https://bitbucket.org/grackle/grackle/raw/default/input/CloudyData_noUVB.h5
tar xzvf LOW.tar.gz
mv LOW/*.dat ./
rmdir LOW
rm LOW.tar.gz

wget --no-check-certificate -O ./HI.tar.gz https://www.dropbox.com/scl/fi/vld7sq9aq2bs5r0wyjgmg/HI.tar.gz?rlkey=ht3yri3x5157xqrxvw94rli0q&st=dc8ig863&dl=1
wget https://bitbucket.org/grackle/grackle/raw/default/input/CloudyData_UVB=HM2012.h5
tar xzvf HI.tar.gz
mv HI/*.dat ./
rmdir HI
rm HI.tar.gz

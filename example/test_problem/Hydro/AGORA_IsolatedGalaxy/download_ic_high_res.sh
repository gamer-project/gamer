wget --no-check-certificate -O ./HI.tar.gz https://www.dropbox.com/scl/fi/m8a6iutr4uf5p2whjsydw/HI.tar.gz?rlkey=h3gz1z7usxzqhuud1kos4zzff
wget --no-check-certificate  -O ./CloudyData_UVB=HM2012.h5 https://github.com/grackle-project/grackle_data_files/raw/main/input/CloudyData_UVB=HM2012.h5
tar xzvf HI.tar.gz
mv HI/*.dat ./
rm -rf HI
rm HI.tar.gz


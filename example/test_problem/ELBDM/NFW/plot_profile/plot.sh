# $1 : Input    file__Initial   Density
# $2 : Input    file__Final     Density
# $3 : Output   Profile name
cp ../$1 $1
cp ../$2 $2

./GAMER_ExtractProfile  -D -i $1 -S -r 1 -n 128 -x 1.5 -y 1.5 -z 1.5 -p
mv AveParDens AveParDens_init
./GAMER_ExtractProfile  -D -i $2 -S -r 1 -n 128 -x 1.5 -y 1.5 -z 1.5 -p
mv AveParDens AveParDens_final
gnuplot plot_profile.gpt
mv Density_Profile.png $3
#! /bin/bash

declare -i t1 t2 t3 t4 t5 t6


# read command-line options
if [ "$1" = "" ]; then
   echo "Please provide the input filename (without the suffix \"_XXXXXX\") [and the range of data ID] :"
   echo -e "\n   sh $0 InputName [StartID EndID [DeltaID]]\n"
   exit 1

elif [ "$2" != "" ] && [ "$3" = "" ]; then
   echo "Please provide the input filename (without the suffix \"_XXXXXX\") [and the range of data ID] :"
   echo -e "\n   sh $0 InputName [StartID EndID [DeltaID]]\n"
   exit 2

elif [ "$2" != "" ] && [ "$3" != "" ]; then
   StartID=$2
   EndID=$3
   AddSuffix=1

   if [ "$4" = "" ]; then dID=1;
   else                   dID=$4;
   fi

else
   StartID=0
   EndID=0
   AddSuffix=0
   dID=1;
fi


# loop over all targeted files
echo ""

for (( t=$StartID; t<=$EndID; t=t+$dID ))
do

#  set the file suffix
   if [ $AddSuffix = 1 ]; then
      t6=t/100000
      t5=(t%100000)/10000
      t4=(t%10000)/1000
      t3=(t%1000)/100
      t2=(t%100)/10
      t1=(t%10)

      SUFFIX="_$t6$t5$t4$t3$t2$t1"
      FILE_IN=$1$SUFFIX
   else
      SUFFIX=""
      FILE_IN=$1
   fi


#  begin sorting
   echo -e "***** Working on data $FILE_IN ... *****"


#  example for loading the center coords from the last row of the file CenterCoords
#  FileCenter=CenterCoords
#  if [ ! -f $FileCenter ]; then
#      echo "ERROR : the file \"$FileCenter\" does not exists !!"
#      exit 4
#  fi

#  x=$(tail -n 1 CenterCoords |awk '{print $4}');
#  y=$(tail -n 1 CenterCoords |awk '{print $5}');
#  z=$(tail -n 1 CenterCoords |awk '{print $6}');

#  ./GAMER_ExtractProfile -i $FILE_IN -o "$SUFFIX" -S -a 1 -r 0.30 -m 8 -a 1 -S -p -R 0.15 -c -x $x -y $y -z $z


#  example for setting all parameters manually
   ./GAMER_ExtractProfile -i $FILE_IN -o "$SUFFIX" -S -a 1 -m 8 -p -r 0.10 -x 0.54 -y 0.64 -z 1.40


   ErrorCode=$?
   if [ $ErrorCode != 0 ]; then
      echo "$0 failed ... (error code = $ErrorCode)"
      exit 3
   fi

   echo -e "***** Working on data $FILE_IN ... done *****\n"

done


exit 0


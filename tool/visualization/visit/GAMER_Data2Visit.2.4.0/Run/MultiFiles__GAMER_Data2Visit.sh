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
   WithSuffix=1

   if [ "$4" = "" ]; then dID=1;
   else                   dID=$4;
   fi

else
   StartID=0
   EndID=0
   WithSuffix=0
   dID=1;
fi


# loop over all targeted files
echo ""

for (( t=$StartID; t<=$EndID; t=t+$dID ))
do

#  set up the file suffix
   if [ $WithSuffix = 1 ]; then
      t6=t/100000
      t5=(t%100000)/10000
      t4=(t%10000)/1000
      t3=(t%1000)/100
      t2=(t%100)/10
      t1=(t%10)

      FILE_IN="$1_$t6$t5$t4$t3$t2$t1"
   else
      FILE_IN="$1"
   fi

   DIR=`echo $FILE_IN | sed 's#.*/##'`


#  begin sorting
   echo -e "***** Working on data $FILE_IN ... *****"

   ./GAMER_Data2Visit -i $FILE_IN -x 0.0 -y 0.0 -z 0.0 -X 1.0 -Y 1.0 -Z 1.0 -m 3 -p -w

   ErrorCode=$?
   if [ $ErrorCode != 0 ]; then
      echo "$0 failed ... (error code = $ErrorCode)"
      exit 3
   fi


#  move data to a new directory
   mkdir $DIR
   mv Level*    $DIR/
   mv Root.silo $DIR/
   mv Box.silo  $DIR/


#  construct the group list for animation
   echo "$DIR/Root.silo" >> Group.visit


   echo -e "***** Working on data $FILE_IN ... done *****\n"

done


exit 0


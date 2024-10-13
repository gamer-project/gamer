#! /bin/bash

declare -i t1 t2 t3 t4 t5 t6


# read command-line options
if [ "$1" = "" ]; then
   echo "Please provide the input filename (without the suffix \"_XXXXXX\") [and the range of data ID] :"
   echo -e "\n   sh $0 InputName [StartID EndID]\n"
   exit 1
elif [ "$2" != "" ] && [ "$3" = "" ]; then
   echo "Please provide the input filename (without the suffix \"_XXXXXX\") [and the range of data ID] :"
   echo -e "\n   sh $0 InputName [StartID EndID]\n"
   exit 2
elif [ "$2" != "" ] && [ "$3" != "" ]; then
   StartID=$2
   EndID=$3
   WithSuffix=1
else
   StartID=0
   EndID=0
   WithSuffix=0
fi


# export variables for gnuplot
export FILE_TEMP=GAMER_Uniform2D2gnuplot__TempFile
export GPT_TERM=png


# loop over all targeted files
echo ""

for (( t=$StartID; t<=$EndID; t++ ))
do

#  set up the file suffix
   if [ $WithSuffix = 1 ]; then
      t6=t/100000
      t5=(t%100000)/10000
      t4=(t%10000)/1000
      t3=(t%1000)/100
      t2=(t%100)/10
      t1=(t%10)

      export FILE_IN="$1_$t6$t5$t4$t3$t2$t1"
   else
      export FILE_IN="$1"
   fi

   export FILE_OUT=`echo $FILE_IN | sed 's#.*/##'`."$GPT_TERM"

   echo "***** Plotting the file \"$FILE_IN\" *****"


#  transform the input 2D data to the gnuplot format
   ./GAMER_Uniform2D2gnuplot -i $FILE_IN -o $FILE_TEMP -n 1

   ErrorCode=$?
   if [ $ErrorCode != 0 ]; then
      echo "GAMER_Uniform2D2gnuplot failed ... (error code = $ErrorCode)"
      exit 3
   fi


#  create figures with gnuplot
   gnuplot gnuplot_pm3d.gpt

   ErrorCode=$?
   if [ $ErrorCode != 0 ]; then
      echo "gnuplot failed ... (error code = $ErrorCode)"
      exit 4
   fi


   echo -e "***** Plot is complete *****\n"

   rm -f $FILE_TEMP

done


unset FILE_TEMP
unset FILE_IN
unset FILE_OUT
unset GPT_TERM

#exit 0


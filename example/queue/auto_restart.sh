#!/bin/bash

declare -i RESTART_COUNT=0 RESTART_COUNT_MAX


RESTART_AUTO=1          # 0/1 <==> auto restart on/off
RESTART_COUNT_MAX=4     # maximum number of auto restarts from the same snapshot


while true; do


#  run gamer
#  ------------------------------------------------------
   mpirun -map-by ppr:1:socket:pe=10 ./gamer
#  ------------------------------------------------------


#  check if the run completes successfully
   ErrorCode=$?
   if [ ${ErrorCode} -eq 0 ]; then
      printf "\n"
      printf "=====================================================================\n"
      printf "   program completed successfully\n"
      printf "=====================================================================\n"
      exit 0

   else
      printf "\n"
      printf "=====================================================================\n"
      printf "   program failed (error code %d) !!\n" $ErrorCode
      printf "=====================================================================\n"

#     restart simulation from the last snapshot automatically
      if [ ${RESTART_AUTO} -eq 1 ]; then

#        check the restart counter
         if [ ${RESTART_COUNT} -lt ${RESTART_COUNT_MAX} ]; then

#           turn on the restart mode in the runtime parameter file
            RESTART_LINE=$(sed -n '/OPT__INIT\ /=' Input__Parameter)
            sed -i "${RESTART_LINE}s/1\ /2\ /" Input__Parameter
            sed -i "${RESTART_LINE}s/3\ /2\ /" Input__Parameter

#           set the restart file
            RESTART_FILE_OLD=$RESTART_FILE_NEW
            RESTART_FILE_NEW=$(ls Data* |tail -n 1)
            ln -sf $RESTART_FILE_NEW RESTART

#           record the restart counter
            if [ "$RESTART_FILE_NEW" = "$RESTART_FILE_OLD" ]; then
               (( RESTART_COUNT ++ ))
            else
               (( RESTART_COUNT = 1 ))
            fi

            printf "\n"
            printf "=====================================================================\n"
            printf "   auto restart (counter %d) ...\n" $RESTART_COUNT
            printf "=====================================================================\n"

         else
            printf "\n"
            printf "=====================================================================\n"
            printf "   exceed the maximum number of auto restart [%d] ... abort !!\n" $RESTART_COUNT_MAX
            printf "=====================================================================\n"
            exit 2
         fi # if [ ${RESTART_COUNT} -lt ${RESTART_COUNT_MAX} ] ... else ...

      else
         exit 1
      fi # if [ ${RESTART_AUTO} -eq 1 ] ... else ...
   fi # if [ ${ErrorCode} -eq 0 ] ... else ...

done # while true

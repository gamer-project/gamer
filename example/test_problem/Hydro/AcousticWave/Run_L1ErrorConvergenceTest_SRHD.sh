#!/bin/bash

PYTHON=python3
TEMP=(HighT LowT)
RES=(032 064 128 256 512)
PREFIX=AcousticWave

sed -i "s/Acoustic_Temp_Bg\ \ \ \ \ 1.0\ \ \ \ /Acoustic_Temp_Bg\ \ \ \ \ 1.0e+10/g" Input__TestProb

for temp in "${TEMP[@]}"; do

   # Create directory
   mkdir ${temp}

   # Change parameters
   if [ ${temp} == "HighT" ]; then
      sed -i "s/Acoustic_Temp_Bg\ \ \ \ \ 1.0e-10/Acoustic_Temp_Bg\ \ \ \ \ 1.0e+10/g" Input__TestProb
   fi
   if [ ${temp} == "LowT" ]; then
      sed -i "s/Acoustic_Temp_Bg\ \ \ \ \ 1.0e+10/Acoustic_Temp_Bg\ \ \ \ \ 1.0e-10/g" Input__TestProb
   fi

   for res in "${RES[@]}"; do
      echo ${temp}_N=${res}

      # Change parameters
      sed -i "s/NX0_TOT_X\ \+[0-9]*\ \+#\ nu/NX0_TOT_X\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ${res}\ \ \ \ \ \ \ \ \ #\ nu/g" Input__Parameter
      sed -i "s/NX0_TOT_Y\ \+[0-9]*\ \+#\ nu/NX0_TOT_Y\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ${res}\ \ \ \ \ \ \ \ \ #\ nu/g" Input__Parameter
      sed -i "s/NX0_TOT_Z\ \+[0-9]*\ \+#\ nu/NX0_TOT_Z\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ${res}\ \ \ \ \ \ \ \ \ #\ nu/g" Input__Parameter

      # Run gamer
      sh clean.sh
      ./gamer 1>>log 2>&1

      # Store the results
      DIR=${temp}/res_${res}
      mkdir ${DIR}
      cp Input__* ${DIR}
      mv log ${DIR}
      mv Record__* ${DIR}
      mv ${PREFIX}_* ${DIR}

   done
done

# Plot the L1 error convergence
${PYTHON} plot_L1error_SRHD.py

echo "Done."

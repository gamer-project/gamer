ENGY_BINS_TABLE="energyBinsTable"
CRindexs=(2.2 2.4 2.6)


# Extract the data from ENGY_BINS_TABLE to the temporary file
temp_file=$(mktemp)
tail -n +2 "$ENGY_BINS_TABLE" > "$temp_file"

for CRindex in "${CRindexs[@]}"
do
    while IFS=$' \t' read -ra line; do
        Emin=${line[0]}
        ObservedEngyEv=$(printf "%9.4e" ${line[1]})
        Emax=${line[2]}
        FileName="Projected_FRB_Data_000035_${ObservedEngyEv}_1.1e6_${CRindex}.h5"

        if [ -f "$FileName" ]; then
            echo "${FileName} exist !!"
        else
            echo "Running with (ObservedEngyEv=${ObservedEngyEv}, CRindex=${CRindex})"
            ./Project FRB_Data_000035.h5 $ObservedEngyEv 1.1e6 $CRindex robitaille_DL07_PAHISMMix.dat >& log-$CRindex-$ObservedEngyEv
        fi
    done < "${temp_file}"
done

rm "${temp_file}"

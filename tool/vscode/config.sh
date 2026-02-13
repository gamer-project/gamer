bin_working=$(cat bin_working)
if [ -z "${bin_working}" ]; then
    echo "Error: Please set the working directory by the task 'set-working-bin' first."
    exit 1
fi

cd ../src
cp ../bin/${bin_working}/generate_make.sh .
sh generate_make.sh
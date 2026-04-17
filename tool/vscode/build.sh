PYTHON=python3

bin_working=$(cat bin_working)
if [ -z "${bin_working}" ]; then
    echo "Error: Please set the working directory by the task 'set-working-bin' first."
    exit 1
fi

cd ../src
make clean

cd ../.vscode
${PYTHON} extract_macros.py

cd ../src
make -j 8

cd ../bin/${bin_working}/
cp ../gamer .

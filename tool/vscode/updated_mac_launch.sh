#!/bin/bash

# Check OS
if [ "$(uname)" == "Darwin" ]; then
    echo "OS: macOS. Start updating launch.json..."
else
    echo "OS: Linux. Do nothing."
    exit 0
fi

cd "$(dirname "$0")"

# If $1 is given, use it as $lldb_mi_path, else prompt the user
if [ -n "$1" ]; then
    lldb_mi_path="$1"
else
    echo -n "Enter your lldb-mi path. (e.g. /Users/user_namer/lldb-mi/build/src/lldb-mi): "
    read lldb_mi_path
fi

# Input file
INPUT_FILE="launch.json"

# Replace MIMode from "gdb" to "lldb"
sed -i '' 's/"MIMode": "gdb"/"MIMode": "lldb"/' "$INPUT_FILE"

# Replace miDebuggerPath to the new path
sed -i '' 's|"miDebuggerPath": "/usr/bin/gdb"|"miDebuggerPath": "$lldb_mi_path"|' "$INPUT_FILE"

# Replace miDebuggerArgs to "-q"
sed -i '' 's/"miDebuggerArgs": "-quiet"/"miDebuggerArgs": "-q"/' "$INPUT_FILE"

echo "launch.json updated successfully."
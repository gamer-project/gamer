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

# Parameters to change
MIMODE="lldb"
MIDEBUGGERPATH=$lldb_mi_path
MIDEBUGGERARGS="-q"

# Input file
INPUT_FILE="launch.json"

# Replace MIMode from "gdb" to "lldb"
if grep -q "\"MIMode\": \"$MIMODE\"" "$INPUT_FILE"; then
    echo "Warning: 'MIMode' is already set to '$MIMODE'."
else
    sed -i '' 's/"MIMode": "gdb"/"MIMode": "'"$MIMODE"'"/' "$INPUT_FILE"
    echo "'MIMode' changed to '$MIMODE'."
fi

# Replace miDebuggerPath to the new path
sed -i '' 's|"miDebuggerPath": .*|"miDebuggerPath": "'"$MIDEBUGGERPATH"'",|' "$INPUT_FILE"
echo "'miDebuggerPath' changed to '$MIDEBUGGERPATH'."


# Replace miDebuggerArgs to "-q"
if grep -q "\"miDebuggerArgs\": \"$MIDEBUGGERARGS\"" "$INPUT_FILE"; then
    echo "Warning: 'miDebuggerArgs' is already set to '$MIDEBUGGERARGS'."
else
    sed -i '' 's/"miDebuggerArgs": "-quiet"/"miDebuggerArgs": "'"$MIDEBUGGERARGS"'"/' "$INPUT_FILE"
    echo "'miDebuggerArgs' changed to '$MIDEBUGGERARGS'."
fi

echo "launch.json updated successfully."
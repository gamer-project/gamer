#!/bin/bash

cd "$(dirname "$0")"

show_help() {
    echo "Usage: $0 [--local | --global] --machine=<name>"
    echo ""
    echo "Options:"
    echo "  --local           Use local settings"
    echo "  --global          Use global settings"
    echo "  --machine=<name>  Specify the machine name"
    echo "  -h, --help        Show this help message"
}

# Parse arguments
LOCAL=false
GLOBAL=false
MACHINE=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --local) LOCAL=true ;;
        --global) GLOBAL=true ;;
        --machine=*) MACHINE="${1#*=}" ;;
        -h|--help) show_help; exit 0 ;;
        --machine) echo "Error: Use --machine=<name> to specify the machine."; exit 1 ;;
        *) echo "Unknown parameter passed: $1"; echo "Use -h or --help for usage information."; exit 1 ;;
    esac
    shift
done

if [ -z "$MACHINE" ]; then
    echo "Error: --machine option is required."
    exit 1
fi

if [ "$LOCAL" = true ] && [ "$GLOBAL" = true ]; then
    echo "Error: Cannot specify both --local and --global."
    exit 1
elif [ "$LOCAL" = true ]; then
    SETTING_FILE="../src/.local_settings"
    SETTING_TYPE="local"
elif [ "$GLOBAL" = true ]; then
    SETTING_FILE="$HOME/.config/gamer/global_settings"
    SETTING_TYPE="global"
else
    echo "Error: Specify either --local or --global."
    exit 1
fi

mkdir -p "$(dirname "$SETTING_FILE")"

# Write the machine name to the setting file
echo "# GAMER setting file" > "$SETTING_FILE"
echo "machine_name          $MACHINE" >> "$SETTING_FILE"

# Check if the write was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to write to $SETTING_FILE."
    exit 1
fi

echo "Successfully wrote $SETTING_TYPE settings to $SETTING_FILE for machine='$MACHINE'."

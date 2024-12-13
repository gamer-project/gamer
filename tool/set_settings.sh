#!/bin/bash

cd "$(dirname "$0")"

# Define valid keys
VALID_KEYS=("machine" "user" "organization")

show_help() {
    echo "Usage: $0 [--local | --global] [--key=value ...] [--delete key ...]"
    echo ""
    echo "Options:"
    echo "  --local              Use local settings"
    echo "  --global             Use global settings"
    echo "  --key=value          Set a parameter"
    echo "  --delete, -d key     Delete a parameter"
    echo "  -h, --help           Show this help message"
    echo ""
    echo "Valid keys and their functionalities:"
    echo "  machine        Specify the machine name"
    echo "  user           Specify the user name"
    echo "  organization   Specify the organization name"
}

# Parse arguments
LOCAL=false
GLOBAL=false
DELETE=false
SET=false
declare -A SETTINGS
DELETE_KEYS=()

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --local) LOCAL=true; shift ;;
        --global) GLOBAL=true; shift ;;
        --delete|-d)
            DELETE=true
            shift # Shift past '--delete' or '-d'
            if [[ "$#" -eq 0 ]]; then
                echo "Error: --delete requires at least one key."
                exit 1
            fi
            while [[ "$#" -gt 0 && "${1:0:2}" != "--" ]]; do
                if [[ ! " ${VALID_KEYS[@]} " =~ " $1 " ]]; then
                    echo "Error: Invalid key '$1' for deletion."
                    exit 1
                fi
                DELETE_KEYS+=("$1")
                shift
            done
            ;;
        -h|--help) show_help; exit 0 ;;
        --*)
            SET=true
            PARAM="${1#--}"
            if [[ "$PARAM" != *=* ]]; then
                echo "Error: Invalid format. Use --key=value."
                exit 1
            fi
            KEY="${PARAM%%=*}"
            VALUE="${PARAM#*=}"
            if [[ ! " ${VALID_KEYS[@]} " =~ " $KEY " ]]; then
                echo "Error: Invalid key '$KEY'."
                exit 1
            fi
            SETTINGS["$KEY"]="$VALUE"
            shift
            ;;
        *) echo "Unknown parameter passed: $1"; echo "Use -h or --help for usage information."; exit 1 ;;
    esac
done

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

if [ "$SET" = true ] && [ "$DELETE" = true ]; then
    echo "Error: Cannot set and delete parameters at the same time."
    exit 1
elif [ "$SET" = false ] && [ "$DELETE" = false ]; then
    echo "Error: No parameters to set or delete."
    exit 1
fi

mkdir -p "$(dirname "$SETTING_FILE")"

# Load existing settings
declare -A EXISTING_SETTINGS
if [ -f "$SETTING_FILE" ]; then
    while read -r LINE; do
        [[ "$LINE" =~ ^#.*$ ]] && continue
        KEY="$(echo "$LINE" | awk '{print $1}')"
        VALUE="$(echo "$LINE" | awk '{print $2}')"
        EXISTING_SETTINGS["$KEY"]="$VALUE"
    done < "$SETTING_FILE"
fi

# Update settings
if [ "$SET" = true ]; then
    for KEY in "${!SETTINGS[@]}"; do
        EXISTING_SETTINGS["$KEY"]="${SETTINGS[$KEY]}"
    done
fi

# Delete specified keys
if [ "$DELETE" = true ]; then
    for KEY in "${DELETE_KEYS[@]}"; do
        unset EXISTING_SETTINGS["$KEY"]
    done
fi

# Write updated settings to file
{
    echo "# GAMER setting file"
    for KEY in "${!EXISTING_SETTINGS[@]}"; do
        echo "$KEY          ${EXISTING_SETTINGS[$KEY]}"
    done
} > "$SETTING_FILE"

# Check if the write was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to write to $SETTING_FILE."
    exit 1
fi

echo "Successfully updated $SETTING_TYPE settings in $SETTING_FILE."

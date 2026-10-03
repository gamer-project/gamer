#!/bin/bash

#################### UTILITIES ####################

cd "$(dirname "$0")"
source ./settings.sh
# Variables imported from settings.sh:
# KEY_DESCRIPTIONS
# VALID_KEYS
# MAX_KEY_LENGTH

# Functions imported from settings.sh:
# print_key()
# show_valid_keys()

show_help() {
    echo "Usage:"
    echo "  $0 [--local | --global] <key>"
    echo "  $0 (-h | --help)"
    echo ""
    echo "Options:"
    echo "  --local              Use local settings"
    echo "  --global             Use global settings"
    echo "  -h, --help           Show this help message"
    echo ""
    show_valid_keys
}

##################### PARSER ######################

LOCAL=false
GLOBAL=false
KEY=""

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --local)
            LOCAL=true
            shift ;;
        --global)
            GLOBAL=true
            shift ;;
        -h|--help)
            show_help
            exit 0 ;;
        -*)
            echo "Error: Unknown option '$1'" >&2
            show_help
            exit 1 ;;
        *)
            if [[ -n "$KEY" ]]; then
                echo "Error: Only one key can be specified." >&2
                exit 1
            fi
            KEY="$1"
            shift ;;
    esac
done

#################### VALIDATION ####################

if [[ -z "$KEY" ]]; then
    echo "Error: You must specify a key to retrieve." >&2
    show_help
    exit 1
fi

if [[ ! " ${VALID_KEYS[@]} " =~ " $KEY " ]]; then
    echo "Error: Invalid key '$KEY'." >&2
    show_valid_keys
    exit 1
fi

################### GET SETTINGS ###################

# Helper to extract a setting from output
extract_value() {
    local key="$1"
    grep -E "^$key[[:space:]]+" | awk '{print $2}'
}

FOUND=false

if [ "$LOCAL" = true ] || [ "$GLOBAL" = false ]; then
    local_output=$(sh ./set_settings.sh --local --list)
    value=$(echo "$local_output" | extract_value "$KEY")
    if [[ -n "$value" ]]; then
        echo "$value"
        FOUND=true
    fi
fi

if [ "$FOUND" = false ] && ([ "$GLOBAL" = true ] || [ "$LOCAL" = false ]); then
    global_output=$(sh ./set_settings.sh --global --list)
    value=$(echo "$global_output" | extract_value "$KEY")
    if [[ -n "$value" ]]; then
        echo "$value"
        FOUND=true
    fi
fi

if [ "$FOUND" = false ]; then
    echo "Key '$KEY' is not set in the requested scope(s)." >&2
    exit 1
fi

#!/bin/bash

################### DEFINE KEYS ###################
# The keys can not contain spaces or start with `-`
declare -A KEY_DESCRIPTIONS
KEY_DESCRIPTIONS=(
    ["machine"]="Specify the machine name"
)
###################################################

# List of valid keys
VALID_KEYS=("${!KEY_DESCRIPTIONS[@]}")

# Keys padded with trailing spaces for formatting output
MAX_KEY_LENGTH=0
for KEY in "${VALID_KEYS[@]}"; do
    if [ ${#KEY} -gt $MAX_KEY_LENGTH ]; then
        MAX_KEY_LENGTH=${#KEY}
    fi
done
declare -A PADDED_KEYS
for KEY in "${VALID_KEYS[@]}"; do
    PADDED_KEYS[$KEY]=$(printf "%-${MAX_KEY_LENGTH}s" "$KEY")
done

show_valid_keys() {
    echo "Valid keys and their functionalities:"
    for KEY in "${!KEY_DESCRIPTIONS[@]}"; do
        echo "  ${PADDED_KEYS[$KEY]}  ${KEY_DESCRIPTIONS[$KEY]}"
    done
}

show_help() {
    echo "Usage:"
    echo "  $0 (--local | --global) [-lv] [--<key>=<value> ...]"
    echo "  $0 (--local | --global) [-lvd] [--delete <key> ...]"
    echo "  $0 [-h | --help]"
    echo ""
    echo "Options:"
    echo "  --local              Use local settings"
    echo "  --global             Use global settings"
    echo "  -l, --list           List current settings"
    echo "  -v, --verbose        Show verbose output"
    echo "  --<key>=<value>      Set a parameter"
    echo "  -d, --delete key     Delete a parameter"
    echo "  --clear-all          Clear all parameters"
    echo "  -h, --help           Show this help message"
    echo ""
    show_valid_keys
}

# Parse arguments
LOCAL=false
GLOBAL=false
LIST=false
VERBOSE=false
DELETE=false
SET=false
declare -A SETTINGS
DELETE_KEYS=()

if [ "$#" -eq 0 ]; then
    printf "$(show_help)\n" >&2
    exit 1
fi

parse_set_parameter() {
    local PARAM="$1"
    local KEY="${PARAM%%=*}"
    if [[ ! " ${VALID_KEYS[@]} " =~ " $KEY " ]]; then
        echo "Error: Invalid key '$KEY'." >&2
        printf "$(show_valid_keys)\n" >&2
        exit 1
    fi
    if [[ "$PARAM" != *=* ]]; then
        echo "Error: Invalid format for key '$KEY'. Use --$KEY=<value>." >&2
        exit 1
    fi
    local VALUE="${PARAM#*=}"
    if [[ -z "$VALUE" ]]; then
        echo "Error: Value for key '$KEY' cannot be empty. Use --$KEY=<value>." >&2
        echo "If you want to delete '$KEY', use --delete $KEY." >&2
        exit 1
    fi
    if [[ "$VALUE" =~ \  ]]; then
        echo "Error: Value for key '$KEY' cannot contain spaces." >&2
        exit 1
    fi
    SETTINGS["$KEY"]="$VALUE"
}

parse_delete_keys() {
    keys_processed=0
    if [[ "$#" -eq 0 ]]; then
        echo "Error: -d or --delete requires at least one key." >&2
        exit 1
    fi
    while [[ "$#" -gt 0 && "${1:0:1}" != "-" ]]; do
        if [[ ! " ${VALID_KEYS[@]} " =~ " $1 " ]]; then
            echo "Error: Invalid key '$1' for deletion." >&2
            printf "$(show_valid_keys)\n" >&2
            exit 1
        fi
        DELETE_KEYS+=("$1")
        shift
        keys_processed=$((keys_processed + 1))
    done
}

while [[ "$#" -gt 0 ]]; do
    if [[ "${1:0:1}" = "-" && "${1:0:2}" != "--" ]]; then # Short options
        OPT="${1:1}"
        shift
        while [[ -n "$OPT" ]]; do # Possibly combined
            case "${OPT:0:1}" in
                l) LIST=true ;;
                v) VERBOSE=true ;;
                h) show_help; exit 0 ;;
                d)
                    DELETE=true
                    # Put the remaining combined short options back to the argument list
                    # Warning: This method will only apply if all other options are not order sensitive
                    if [[ -n "${OPT:1}" ]]; then
                        set -- "$@" "-$(echo "${OPT:1}" | tr -d 'd')"
                    fi
                    parse_delete_keys "$@"
                    shift $keys_processed
                    break ;;
                *)
                    echo "Error: Unknown option: -${OPT:0:1}" >&2
                    echo "Use -h or --help for usage information." >&2
                    exit 1 ;;
            esac
            OPT="${OPT:1}"
        done
    else                                       # Long options
        case $1 in
            --local) LOCAL=true; shift ;;
            --global) GLOBAL=true; shift ;;
            --list) LIST=true; shift ;;
            --verbose) VERBOSE=true; shift ;;
            --help) show_help; exit 0 ;;
            --delete)
                DELETE=true
                shift
                parse_delete_keys "$@"
                shift $keys_processed ;;
            --clear-all)
                DELETE=true
                DELETE_KEYS=("${VALID_KEYS[@]}")   # Set DELETE_KEYS to all valid keys
                shift ;;
            --*)
                SET=true
                parse_set_parameter "${1#--}"
                shift ;;
            *)
                echo "Error: Unknown parameter passed: $1" >&2
                echo "Use -h or --help for usage information." >&2
                exit 1 ;;
        esac
    fi
done

cd "$(dirname "$0")"

# Ensure either --local or --global is specified
if [ "$LOCAL" = true ] && [ "$GLOBAL" = true ]; then
    echo "Error: Cannot specify both --local and --global." >&2
    exit 1
elif [ "$LOCAL" = true ]; then
    SETTING_FILE="../src/.local_settings"
    SETTING_TYPE="local"
elif [ "$GLOBAL" = true ]; then
    SETTING_FILE="$HOME/.config/gamer/global_settings"
    SETTING_TYPE="global"
else
    echo "Error: Specify either --local or --global." >&2
    exit 1
fi

# Ensure mutually exclusive operations
if [ "$SET" = true ] && [ "$DELETE" = true ]; then
    echo "Error: Cannot set and delete parameters at the same time." >&2
    exit 1
elif [ "$SET" = false ] && [ "$DELETE" = false ]; then
    LIST=true
fi

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

# Create the directory and the settings file if it doesn't exist
DIR="$(dirname "$SETTING_FILE")"
if [ ! -d "$DIR" ]; then
    mkdir -p "$DIR"
    [ "$VERBOSE" = true ] && echo "Created directory '$(realpath "$DIR")'."
fi
if [ ! -f "$SETTING_FILE" ]; then
    if touch "$SETTING_FILE"; then
        [ "$VERBOSE" = true ] && echo "Created file $(basename "$SETTING_FILE") in '$(realpath "$DIR")'."
    else
        echo "Fatal: Failed to create file $(basename "$SETTING_FILE") in '$(realpath "$DIR")'." >&2
        exit 1
    fi
fi

# The head of message to the user
if [ "$LIST" = true ]; then # Header for listing settings
    echo "$SETTING_TYPE settings in $(realpath "$SETTING_FILE")"
    echo "---------------------------"
elif [ "$VERBOSE" = true ]; then
    echo "Loaded $SETTING_TYPE settings from '$(realpath "$SETTING_FILE")'."
fi

# Main loop to update or delete settings
VERBOSE_OUTPUT="" # Store the echo output for verbose mode
for KEY in "${VALID_KEYS[@]}"; do
    OLD_VALUE="${EXISTING_SETTINGS[$KEY]}"
    NEW_VALUE="${SETTINGS[$KEY]}"
    PADDED_KEY="${PADDED_KEYS[$KEY]}"

    if [ "$SET" = true ] && [ -n "$NEW_VALUE" ]; then # The key will be set

        EXISTING_SETTINGS["$KEY"]="$NEW_VALUE" # Update or set new value
        if [ "$LIST" = true ]; then
            if [ -z "$OLD_VALUE" ]; then
                echo "$PADDED_KEY   $NEW_VALUE (new)"
            else
                echo "$PADDED_KEY   $OLD_VALUE -> $NEW_VALUE"
            fi
        fi
        if [ "$VERBOSE" = true ]; then
            if [ -z "$OLD_VALUE" ]; then
                VERBOSE_OUTPUT+="Set '$KEY' to '$NEW_VALUE'.\n"
            else
                VERBOSE_OUTPUT+="Updated '$KEY' from '$OLD_VALUE' to '$NEW_VALUE'.\n"
            fi
        fi

    elif [ "$DELETE" = true ] && [[ " ${DELETE_KEYS[@]} " =~ " $KEY " ]]; then # The key will be deleted

        unset EXISTING_SETTINGS["$KEY"] # Delete the key
        if [ "$LIST" = true ] && [ -n "$OLD_VALUE" ]; then
            echo "$PADDED_KEY   $OLD_VALUE -> (deleted)"
        fi
        if [ "$VERBOSE" = true ] && [ -n "$OLD_VALUE" ]; then
            VERBOSE_OUTPUT+="Deleted '$KEY' which was set to '$OLD_VALUE'.\n"
        fi

    elif [ "$LIST" = true ] && [ -n "$OLD_VALUE" ]; then
        echo "$PADDED_KEY   $OLD_VALUE"
    fi
done
[ "$LIST" = true ] && echo ""

# Write updated settings to file
{
    echo "# GAMER setting file"
    for KEY in "${!EXISTING_SETTINGS[@]}"; do
        echo "${PADDED_KEYS[$KEY]}    ${EXISTING_SETTINGS[$KEY]}"
    done
} > "$SETTING_FILE"

# Check if writing to file was successful
if [ $? -ne 0 ]; then
    echo "Fatal: Failed to write to '$SETTING_FILE'." >&2
    exit 1
elif [ "$VERBOSE" = true ]; then
    echo -en "$VERBOSE_OUTPUT"
fi
if [ "$SET" = true ] || [ "$DELETE" = true ]; then
    echo "Successfully updated $SETTING_TYPE settings."
fi

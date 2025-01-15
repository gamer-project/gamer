#!/bin/bash

################### DEFINE KEYS ###################
# The keys can not contain spaces or start with `-`
declare -A KEY_DESCRIPTIONS
KEY_DESCRIPTIONS=(
    ["machine"]="Specify the machine name"
)
###################################################

#################### UTILITIES ####################

# List of valid keys
VALID_KEYS=("${!KEY_DESCRIPTIONS[@]}")

# For padding keys with trailing spaces to format output
MAX_KEY_LENGTH=0
for key in "${VALID_KEYS[@]}"; do
    if [ ${#key} -gt $MAX_KEY_LENGTH ]; then
        MAX_KEY_LENGTH=${#key}
    fi
done

# Print keys in a formatted way
print_key() {
    # $1 : the key name
    # $2 : the key value or additional message
    # $3 : indent number (optional, default 0)
    printf "%${3}s" ""
    printf "%-${MAX_KEY_LENGTH}s   %s\n" "$1" "$2"
}

show_valid_keys() {
    echo "Valid keys and their functionalities:"
    for key in "${!KEY_DESCRIPTIONS[@]}"; do
        print_key "$key" "${KEY_DESCRIPTIONS[$key]}" 2
    done
}

show_help() {
    echo "Usage:"
    echo "  $0 (--local | --global) [-l] --<key>=<value> ..."
    echo "  $0 (--local | --global) [-l] (-d | --delete) <key> ..."
    echo "  $0 (--local | --global) (-l | --list)"
    echo "  $0 (-h | --help)"
    echo ""
    echo "Options:"
    echo "  --local              Use local settings"
    echo "  --global             Use global settings"
    echo "  --<key>=<value>      Set a parameter"
    echo "  -d, --delete key     Delete a parameter"
    echo "  --clear-all          Clear all parameters"
    echo "  -l, --list           List current settings"
    echo "  -h, --help           Show this help message"
    echo ""
    show_valid_keys
}

##################### PARSER ######################

# Parser parameters
LOCAL=false
GLOBAL=false
LIST=false
SET=false
DELETE=false
KEYS_PROCESSED=0
declare -A SETTINGS
DELETE_KEYS=()

# Parser utility functions
parse_set_parameter() {
    local param="$1"
    local key="${param%%=*}"
    if [[ "$param" != *=* ]]; then
        echo "Error: Invalid format for the key '$key'. Use --$key=<value>." >&2
        exit 1
    fi
    local value="${param#*=}"
    SETTINGS["$key"]="$value"
}

parse_delete_keys() {
    KEYS_PROCESSED=0
    while [[ "$#" -gt 0 && "${1:0:1}" != "-" ]]; do
        DELETE_KEYS+=("$1")
        shift
        KEYS_PROCESSED=$((KEYS_PROCESSED + 1))
    done
}

# Main parser loop
while [[ "$#" -gt 0 ]]; do
    if [[ "${1:0:1}" = "-" && "${1:0:2}" != "--" ]]; then # Short options
        opt="${1:1}"
        shift
        while [[ -n "$opt" ]]; do # Possibly combined
            case "${opt:0:1}" in
                l) LIST=true ;;
                h) show_help; exit 0 ;;
                d)
                    DELETE=true
                    # Put the remaining combined short options back to the argument list
                    # Warning: This method will only apply if all other options are not order sensitive
                    if [[ -n "${opt:1}" ]]; then
                        set -- "$@" "-$(echo "${opt:1}" | tr -d 'd')"
                    fi
                    parse_delete_keys "$@"
                    shift $KEYS_PROCESSED
                    break ;;
                *)
                    echo "Error: Unknown option: -${opt:0:1}" >&2
                    printf "$(show_help)\n" >&2
                    exit 1 ;;
            esac
            opt="${opt:1}"
        done
    else                                       # Long options
        case $1 in
            --local)
                LOCAL=true
                SETTING_TYPE="local"
                SETTING_FILE="../../src/.local_settings"
                shift ;;
            --global)
                GLOBAL=true;
                SETTING_TYPE="global"
                SETTING_FILE="$HOME/.config/gamer/global_settings"
                shift ;;
            --list) LIST=true; shift ;;
            --help) show_help; exit 0 ;;
            --delete)
                DELETE=true
                shift
                parse_delete_keys "$@"
                shift $KEYS_PROCESSED ;;
            --clear-all)
                DELETE=true
                DELETE_KEYS=("${VALID_KEYS[@]}")   # Set DELETE_KEYS to all valid keys
                shift ;;
            --*)
                SET=true
                parse_set_parameter "${1#--}"
                shift ;;
            *)
                echo "Error: Unknown option: $1" >&2
                printf "$(show_help)\n" >&2
                exit 1 ;;
        esac
    fi
done

############### VALIDATE PARAMETERS ###############

# Ensure at least one operation is specified
if [ "$SET" = false ] && [ "$DELETE" = false ] && [ "$LIST" = false ]; then
    echo "Error: Specify at least one operation." >&2
    printf "$(show_help)\n" >&2
    exit 1
fi

# Validate the keys and values for setting
if [ "$SET" = true ]; then
    for key in "${!SETTINGS[@]}"; do
        if [[ ! " ${VALID_KEYS[@]} " =~ " $key " ]]; then
            echo "Error: Invalid key '$key'." >&2
            printf "$(show_valid_keys)\n" >&2
            exit 1
        fi
        if [[ -z "${SETTINGS[$key]}" ]]; then
            echo "Error: The value for the key '$key' cannot be empty. Use --$key=<value>." >&2
            exit 1
        fi
        if [[ "${SETTINGS[$key]}" =~ \  ]]; then
            echo "Error: The value for the key '$key' cannot contain spaces." >&2
            exit 1
        fi
    done
fi

# Ensure mutually exclusive operations
if [ "$SET" = true ] && [ "$DELETE" = true ]; then
    echo "Error: Cannot set and delete parameters at the same time." >&2
    exit 1
fi

# Validate the keys for deletion
if [ "$DELETE" = true ]; then
    if [ ${#DELETE_KEYS[@]} -eq 0 ]; then
        echo "Error: No keys specified for deletion." >&2
        exit 1
    fi
    for key in "${DELETE_KEYS[@]}"; do
        if [[ ! " ${VALID_KEYS[@]} " =~ " $key " ]]; then
            echo "Error: Invalid key '$key' for deletion." >&2
            printf "$(show_valid_keys)\n" >&2
            exit 1
        fi
    done
fi

# Ensure either --local or --global is specified
if [ "$LOCAL" = true ] && [ "$GLOBAL" = true ]; then
    echo "Error: Cannot specify both --local and --global." >&2
    exit 1
elif [ "$LOCAL" = false ] && [ "$GLOBAL" = false ]; then
    echo "Error: Specify either --local or --global." >&2
    exit 1
fi

################ LOAD SETTINGS FILE ###############

cd "$(dirname "$0")"

# Load if the settings file exists
declare -A EXISTING_SETTINGS
if [ -f "$SETTING_FILE" ]; then
    while read -r LINE; do
        [[ "$LINE" =~ ^#.*$ ]] && continue
        key="$(echo "$LINE" | awk '{print $1}')"
        value="$(echo "$LINE" | awk '{print $2}')"
        EXISTING_SETTINGS["$key"]="$value"
    done < "$SETTING_FILE"
fi

################# UPDATE SETTINGS #################

# The head of the list
if [ "$LIST" = true ]; then # Header for listing settings
    echo "$SETTING_TYPE settings in $SETTING_FILE"
    echo "---------------------------"
fi

# Main loop to update or delete settings
for key in "${VALID_KEYS[@]}"; do
    old_value="${EXISTING_SETTINGS[$key]}"
    new_value="${SETTINGS[$key]}"

    if [ "$SET" = true ] && [ -n "$new_value" ]; then # The key will be set

        EXISTING_SETTINGS["$key"]="$new_value" # Update or set new value
        if [ "$LIST" = true ]; then
            if [ -z "$old_value" ]; then
                print_key "$key" "$new_value (new)"
            else
                print_key "$key" "$old_value -> $new_value"
            fi
        fi

    elif [ "$DELETE" = true ] && [[ " ${DELETE_KEYS[@]} " =~ " $key " ]]; then # The key will be deleted

        unset EXISTING_SETTINGS["$key"] # Delete the key
        if [ "$LIST" = true ] && [ -n "$old_value" ]; then
            print_key "$key" "$old_value -> (deleted)"
        fi

    elif [ "$LIST" = true ] && [ -n "$old_value" ]; then # List the existing settings
        print_key "$key" "$old_value"
    fi
done
[ "$LIST" = true ] && echo "---------------------------"

################ SAVE SETTINGS FILE ###############

if [ "$SET" = true ] || [ "$DELETE" = true ]; then

    # Create the directory and the settings file if it doesn't exist
    dir="$(dirname "$SETTING_FILE")"
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Created directory $dir."
    fi
    if [ ! -f "$SETTING_FILE" ]; then
        if touch "$SETTING_FILE"; then
            echo "Created file $(basename "$SETTING_FILE") in $dir."
        else
            echo "Fatal: Failed to create file $(basename "$SETTING_FILE") in $dir." >&2
            exit 1
        fi
    fi

    # Write updated settings to file
    {
        echo "# GAMER setting file"
        for key in "${!EXISTING_SETTINGS[@]}"; do
            print_key "${key}" "${EXISTING_SETTINGS[$key]}"
        done
    } > "$SETTING_FILE"

    # Check if writing to file was successful
    if [ $? -ne 0 ]; then
        echo "Fatal: Failed to write to '$SETTING_FILE'." >&2
        exit 1
    fi
    echo "Successfully updated $SETTING_TYPE settings."
fi

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

###################################################

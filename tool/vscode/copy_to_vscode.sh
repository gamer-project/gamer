#!/bin/bash

# Change to the directory where the script is located
cd "$(dirname "$0")"

# Set source and destination paths relative to the script's location
SOURCE_DIR="."
TARGET_DIR="../../.vscode"

# Specify files to exclude (space-separated list, e.g., "file1 file2")
EXCLUDE_FILES="copy_to_vscode.sh README.md"

if [ ! -d "$TARGET_DIR" ]; then
  echo "Target directory $TARGET_DIR does not exist."
  mkdir -p "$TARGET_DIR"
  echo "Target directory $TARGET_DIR created."
fi

# Function to check if a file is in the exclude list
is_excluded() {
  for excluded in $EXCLUDE_FILES; do
    if [ "$excluded" = "$1" ]; then
      return 0
    fi
  done
  return 1
}

for file in "$SOURCE_DIR"/*; do
  filename=$(basename "$file")

  if is_excluded "$filename"; then
    continue
  fi

  if [ -e "$TARGET_DIR/$filename" ]; then
    echo "File $filename already exists in $TARGET_DIR. Overwrite? (y/n)"
    read -r answer

    if [ "$answer" = "y" ] || [ "$answer" = "Y" ]; then
      cp "$file" "$TARGET_DIR/$filename"
      echo "$filename overwritten."
    else
      echo "$filename not copied."
    fi
  else
    cp "$file" "$TARGET_DIR/$filename"
    echo "$filename copied."
  fi

done

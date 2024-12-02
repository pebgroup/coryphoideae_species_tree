#!/bin/bash

# Check if the user provided an input folder as the first argument
if [[ -z "$1" ]]; then
  echo "Usage: $0 /path/to/input_folder /path/to/output_folder"
  exit 1
fi

# Use the provided input folder path
input_folder="$1"

# Check if the user provided an output folder as the second argument
if [[ -z "$2" ]]; then
  echo "Usage: $0 /path/to/input_folder /path/to/output_folder"
  exit 1
fi

# Use the provided output folder path
output_folder="$2"

# Ensure the output folder exists
mkdir -p "$output_folder"

# Loop through all .tre files in the input folder
for file in "$input_folder"/*.tre; do
  # Check if the file exists (handles cases with no matching files)
  if [[ -e "$file" ]]; then
    echo "Processing file: $file"
    # Example: Copy the file to the output folder
    nw_reroot -l "$file" 1079 1080 1081 1082 > "$output_folder"
  else
    echo "No .tre files found in $input_folder"
    break
  fi
done

import re
import json

'''
This script updates the "defines" section in the .vscode/c_cpp_properties.json file in Visual Studio Code
by directly parsing the Makefile.log file to extract relevant compiler macros (definitions). It provides an automated
way to synchronize project-specific preprocessor definitions with the VSCode configuration. This approach allows
VSCode to recognize the defines directly from your build configuration, improving VSCode IntelliSense and error
detection without manual intervention.

Usage:
1. Place this script under the `.vscode` folder of your project.
2. Set the paths to `Makefile.log` and `c_cpp_properties.json` in the script as follows:
   - `makefile_log_path`: Path to `Makefile.log`, which contains compiler settings output (default: "../src/Makefile.log").
   - `c_cpp_properties_path`: Path to the VSCode C++ configuration file (default: "c_cpp_properties.json").
3. Run the script each time `Makefile.log` is updated, or set it as a pre-build or post-build task in VSCode
   to keep the configuration in sync.
'''
# Path to Makefile.log and c_cpp_properties.json
makefile_log_path = "../src/Makefile.log"
c_cpp_properties_path = "c_cpp_properties.json"

# Pattern to match the setting in the format of " MODEL : HYDRO"
pattern = re.compile(r":\s+(\w+)\s*:\s+(\w+)")
defines = []

# Read Makefile.log and extract macros
with open(makefile_log_path, 'r') as log_file:
    print(f"Reading {makefile_log_path} and extracting defines...")
    for line in log_file:
        match = pattern.search(line)
        if match:
            key, value = match.groups()
            if value == 'False':
                continue
            elif value == 'True':
                defines.append(f"{key}")
            else:
                defines.append(f"{key}={value}")
print(f"Extracted {len(defines)} macros from {makefile_log_path}.")

# Load c_cpp_properties.json
with open(c_cpp_properties_path, 'r') as cpp_properties_file:
    print(f"Loading {c_cpp_properties_path}...")
    cpp_properties = json.load(cpp_properties_file)

# Update the "defines" array in c_cpp_properties.json
print("Updating defines in c_cpp_properties.json...")
cpp_properties['configurations'][0]['defines'] = defines

# Write the updated content back to c_cpp_properties.json
with open(c_cpp_properties_path, 'w') as cpp_properties_file:
    print(f"Writing updates to {c_cpp_properties_path}...")
    json.dump(cpp_properties, cpp_properties_file, indent=4)

print("Update completed successfully.")

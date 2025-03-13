# Integrating GAMER with Visual Studio Code

Run `sh tool/vscode/copy_to_vscode.sh` to copy necessary files to the `.vscode` directory to integrate GAMER with VS Code.

Please check the [wiki](https://github.com/gamer-project/gamer/wiki/Developing-with-VS-Code) for more information.

## Files in this directory
- `c_cpp_properties.json`: Contains the IntelliSense configuration for VS Code.
- `gamercpp.natvis`: Contains the visualization configuration for VS Code debugger.
- `launch.json`: Contains the debug configuration for VS Code.
- `settings.json`: Contains the settings for the editor in VS Code.
- `tasks.json`: Contains the build configuration for VS Code.
- `extract_macros.py`: Script to extract the macros from the `Makefile.log` to `c_cpp_properties.json`.
- `copy_to_vscode.sh`: Script to copy the above files to `.vscode` directory.
- `bin_working`: File for storing the name of the working directory under `bin/`.
- `set_bin_working.sh`: Script to set up the path of the input files and the executable.
- `build.sh`: Script for the build task.
- `clean_work_dir.sh`: Script for the clean task.
- `config.sh`: Script for the configuration task.
- `README.md`: This file.

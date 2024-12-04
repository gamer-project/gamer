# Integrating Gamer with Visual Studio Code

## Setup
- Install VS Code Extension: "C/C++"
- To disable auto-formatting, add the following to your `settings.json`:
    ```json
    {
        "C_Cpp.formatting": "disabled"
    }
    ```

## Usage
- Run `sh tool/vscode/copy_to_vscode.sh` to copy the necessary files to `.vscode` directory to integrate Gamer with VS Code.
- After configuring Gamer with `configure.py`, press `Ctrl + Shift + B` to build gamer. This will update the macros, providing IntelliSense highlighting support.
- Press `F5` to debug Gamer.

## Files in this directory
- `c_cpp_properties.json`: Contains the IntelliSense configuration for VS Code.
- `gamercpp.natvis`: Contains the visualization configuration for VS Code debugger.
- `launch.json`: Contains the debug configuration for VS Code.
- `tasks.json`: Contains the build configuration for VS Code.
- `extract_macros.py`: Script to extract the macros from the `Makefile.log` to `c_cpp_properties.json`.
- `copy_to_vscode.sh`: Script to copy the above files to `.vscode` directory and set up the working path.
- `README.md`: This file.

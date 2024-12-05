# Developing Gamer with Visual Studio Code

This guide provides step-by-step instructions on how to set up and use Visual Studio Code (VS Code) for developing the Gamer codebase.

## Setup

### Prerequisites

- **Quick Start**: Follow the instructions in the [[Quick Start | Quick-Start]] guide to download Gamer and do at least one [[demo | Quick-Start:-1D-Shock-Tube]].
- **Visual Studio Code**: Download and install from [https://code.visualstudio.com/](https://code.visualstudio.com/).
- **C/C++ Extension**: Install the "C/C++" extension from the VS Code Marketplace.

### Setting up the workspace

1. **Launch VS Code**.
2. **Open the Gamer Project Folder**:
   - Go to `File` > `Open Folder...`.
   - Select your Gamer project directory.

### Disabling auto-formatting

To disable auto-formatting in VS Code, add the following to your `settings.json`:

Access `settings.json` by navigating to `File` > `Preferences` > `Settings`, then search for `C_Cpp.formatting` and set it to `"disabled"`.

### Configuring VS Code to integrate with Gamer

Run the following script from the root directory of the Gamer project:
```bash
sh tool/vscode/copy_to_vscode.sh
```
This script copies the necessary configuration files to the `.vscode` directory, integrating Gamer with VS Code. This script will ask you to enter the working directory under `bin/` where the executables are located. These paths of working directories are located in the `launch.json` and `tasks.json` files.

## Developing with VS Code

### Run build task

After configuring Gamer with [configure.py](https://github.com/gamer-project/gamer/wiki/Installation%3A-Configure.py), select `Terminal` > `Run Build Task...` or press `Ctrl + Shift + B` to start the build task. This will update the macros, providing IntelliSense highlighting support.

Note: To make sure the debugger works correctly, ensure the compiler flags in `Makefile` are set to `-g -O0`. (TBD: Add a argument to `configure.py` to set the flags.)

### Start debugging

To start debugging Gamer, select `Run` > `Start Debugging` or press `F5`. This will launch the debugger. See the [official documentation](https://code.visualstudio.com/docs/editor/debugging) to learn more about debugging with VS Code.

## Understanding Configuration Files

These are the configuration files that are copied to the `.vscode` directory:
- `c_cpp_properties.json`: Configures IntelliSense settings, including include paths and macros. See the [official documentation](https://code.visualstudio.com/docs/cpp/c-cpp-properties-schema-reference) for more information and the [document of IntelliSense](https://code.visualstudio.com/docs/editor/intellisense) to learn about IntelliSense in VS Code.
- `launch.json`: Contains debugging configurations such as executable paths and arguments. See the [official documentation](https://code.visualstudio.com/docs/cpp/launch-json-reference) for more information.
- `tasks.json`: Specifies the build and other dependent tasks. Learn more about [tasks in VS Code](https://code.visualstudio.com/docs/editor/tasks).
- `gamercpp.natvis`: Defines custom visualizations for data structures in the debugger. Learn more about [customizing the visualization of data structures](https://learn.microsoft.com/en-us/visualstudio/debugger/create-custom-views-of-native-objects?view=vs-2022).

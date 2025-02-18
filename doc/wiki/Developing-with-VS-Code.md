# Developing GAMER with Visual Studio Code

This guide provides step-by-step instructions on how to set up and use Visual Studio Code (VS Code) for developing the GAMER codebase.

## Setup

### Prerequisites

- **Quick Start**: Follow the instructions in the [[Quick Start | Quick-Start]] guide to download GAMER and do at least one [[demo | Quick-Start:-1D-Shock-Tube]].
- **Visual Studio Code**: Download and install from [https://code.visualstudio.com/](https://code.visualstudio.com/).
- **C/C++ Extension**: Install the "C/C++" extension from the VS Code Marketplace.

### Setting up the workspace

1. **Launch VS Code**.
2. **Open the GAMER Project Folder**:
   - Go to `File` > `Open Folder...`.
   - Select your GAMER project directory.

### Disabling auto-formatting

To disable auto-formatting in VS Code, add the following to your `settings.json`:

Access `settings.json` by navigating to `File` > `Preferences` > `Settings`, then search for `C_Cpp.formatting` and set it to `"disabled"`.

### Configuring VS Code to integrate with GAMER

Run the following script from the root directory of the GAMER project:
```bash
sh tool/vscode/copy_to_vscode.sh
```
This script copies the necessary configuration files to the `.vscode` directory, integrating GAMER with VS Code.

## Developing with VS Code

Before starting to run the tasks below, **set the working directory** by selecting `Terminal` > `Run Task...` > `set-working-bin` and entering the name of the working directory under `bin/` where the input files are located.

### Configure GAMER

Select `Terminal` > `Run Task...` > `build-GAMER` to configure GAMER with the `generate_make.sh` script in your working directory.

### Build

After configuring GAMER with [configure.py](https://github.com/gamer-project/gamer/wiki/Installation%3A-Configure.py) or `generate_make.sh`, select `Terminal` > `Run Build Task...` or press `Ctrl + Shift + B` to start the build task. This will update the macros, providing IntelliSense highlighting support.

> [!IMPORTANT]  
> To make sure the debugger works correctly, ensure the compiler flags in `Makefile` are set to `-g -O0`. (TBD: Add a argument to `configure.py` to set the flags.)

### Start debugging

To start debugging GAMER, select `Run` > `Start Debugging` or press `F5`. After enter the working directory, the debugger will be launched. See the [official documentation](https://code.visualstudio.com/docs/editor/debugging) to learn more about debugging with VS Code.

> [!NOTE]  
> If `gdb` is not supported on macOS, you can set up `lldb` as the debugger. Make sure [`lldb-mi`](https://github.com/lldb-tools/lldb-mi) is installed. Then select `Terminal` > `Run Task...` > `updated_mac_launch`. This task updates the debugger path in `launch.json` to your `lldb-mi` installation.
> For manual setup or additional details, refer to the [official documentation](https://code.visualstudio.com/docs/cpp/launch-json-reference).

### Clean working directory

Select `Terminal` > `Run Task...` > `clean-work-dir` to clean the working directory with the `clean.sh` script in your working directory.

## Understanding Configuration Files

These are the configuration files that are copied to the `.vscode` directory:
- `c_cpp_properties.json`: Configures IntelliSense settings, including include paths and macros. See the [official documentation](https://code.visualstudio.com/docs/cpp/c-cpp-properties-schema-reference) for more information and the [document of IntelliSense](https://code.visualstudio.com/docs/editor/intellisense) to learn about IntelliSense in VS Code.
- `launch.json`: Contains debugging configurations such as executable paths and arguments. See the [official documentation](https://code.visualstudio.com/docs/cpp/launch-json-reference) for more information.
- `tasks.json`: Specifies the build and other dependent tasks. Learn more about [tasks in VS Code](https://code.visualstudio.com/docs/editor/tasks).
- `gamercpp.natvis`: Defines custom visualizations for data structures in the debugger. Learn more about [customizing the visualization of data structures](https://learn.microsoft.com/en-us/visualstudio/debugger/create-custom-views-of-native-objects?view=vs-2022).

This guide provides step-by-step instructions for setting up and using Visual Studio Code (VS Code) to develop the GAMER codebase.

## Setup

### Prerequisites

- **Quick Start**: Follow the [[Quick Start | Quick-Start]] guide to download GAMER and complete at least one [[demo | Quick-Start:-1D-Shock-Tube]].
- **Visual Studio Code**: Download and install VS Code from [https://code.visualstudio.com/](https://code.visualstudio.com/).
- **C/C++ Extension**: Install the "C/C++" extension from the VS Code Marketplace.

### Setting Up the Workspace

1. **Launch VS Code**.
2. **Open the GAMER Project Folder**:
   - Go to `File` > `Open Folder...`.
   - Select your GAMER project directory.

> [!TIP]
> When using remote-SSH, open the directory as an absolute path to avoid [this issue](https://github.com/microsoft/vscode-cpptools/issues/4818).

### Configuring VS Code for GAMER

Run the following script from the root directory of the GAMER project:
```bash
sh tool/vscode/copy_to_vscode.sh
```
This script copies the necessary configuration files to the `.vscode` directory, integrating GAMER with VS Code.

## Developing with VS Code

Before running the tasks below, **set the working directory** by selecting `Terminal` > `Run Task...` > `set-working-bin` and entering the name of the working directory under `bin/` where the input files are located.

### Configuring GAMER

Select `Terminal` > `Run Task...` > `config-GAMER` to configure GAMER using the `generate_make.sh` script in your working directory.

### Building GAMER

After configuring GAMER with [configure.py](https://github.com/gamer-project/gamer/wiki/Installation%3A-Configure.py) or `generate_make.sh`, select `Terminal` > `Run Task...` > `build-GAMER` to start the build process. This updates the macros and enables IntelliSense highlighting.

> [!TIP]
> To configure and build GAMER in one step, select `Terminal` > `Run Build Task...` or press `Ctrl + Shift + B` to run `config-GAMER` and `build-GAMER` sequentially.

### Debugging GAMER

To start debugging, select `Run` > `Start Debugging` or press `F5`. After entering the working directory, the debugger will launch. See the [official documentation](https://code.visualstudio.com/docs/editor/debugging) to learn more about debugging in VS Code.

> [!IMPORTANT]
> Ensure the compiler flags in `Makefile` are set to `-g -O0` for debugging. (TBD: Add an argument to `configure.py` to set these flags.)

> [!NOTE]
> If `gdb` is not supported on macOS, you can set up `lldb` as the debugger. Ensure [`lldb-mi`](https://github.com/lldb-tools/lldb-mi) is installed, then select `Terminal` > `Run Task...` > `updated_mac_launch`. This updates the debugger path in `launch.json` to your `lldb-mi` installation.
> For manual setup or additional details, refer to the [official documentation](https://code.visualstudio.com/docs/cpp/launch-json-reference).

### Cleaning the Working Directory

Select `Terminal` > `Run Task...` > `clean-work-dir` to clean the working directory using the `clean.sh` script in your working directory.

## Understanding Configuration Files

The following configuration files are copied to the `.vscode` directory:
- `c_cpp_properties.json`: Configures IntelliSense settings, including include paths and macros. See the [schema reference](https://code.visualstudio.com/docs/cpp/c-cpp-properties-schema-reference) and [IntelliSense documentation](https://code.visualstudio.com/docs/editor/intellisense) for details.
- `launch.json`: Defines debugging configurations such as executable paths and arguments. See the [official documentation](https://code.visualstudio.com/docs/cpp/launch-json-reference) for more information.
- `settings.json`: Specifies editor settings, such as indentation spaces and file types for extensions. See the [VS Code settings guide](https://code.visualstudio.com/docs/editor/settings) for details.
- `tasks.json`: Defines build and auxiliary tasks. Learn more about [tasks in VS Code](https://code.visualstudio.com/docs/editor/tasks).
- `gamercpp.natvis`: Customizes data structure visualizations in the debugger. Learn more about [customizing native object views](https://learn.microsoft.com/en-us/visualstudio/debugger/create-custom-views-of-native-objects?view=vs-2022).

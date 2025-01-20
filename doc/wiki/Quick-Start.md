### Download the Code

Follow the instructions in [[Download]] to get the latest version of GAMER.

### Setup Your Machine Configuration File

The machine configuration file is located under `configs`. This file contains the library paths, compiler types, compilation flags, and GPU compute capability. Check the existing configuration files to see if one matches your machine.

If no suitable configuration file is available, you will need to create a new one. Follow the instructions in [[Machine Configuration File | Installation:-Machine-Configuration-File]] to set up a configuration file for your machine.

If a configuration file matches your machine, you can set it as the default by running

```bash
sh tool/config/set_settings.sh --local --machine=your_machine
```

For example, setting `--machine=pleiades` with the above command will use the `configs/pleiades.config` machine configuration when compiling the code.

### Quick Demos

Read the following guides to learn how to configure, compile, and run GAMER.

1. [[1D Shock Tube: OpenMP with/without GPU acceleration | Quick-Start:-1D-Shock-Tube]]
2. [[3D Blast Wave: hybrid MPI/OpenMP/GPU + yt analysis | Quick-Start:-3D-Blast-Wave]]

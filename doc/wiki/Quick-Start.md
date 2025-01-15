### Download the Code

Follow the instructions in [[Download]] to get the latest version of GAMER.

### Setup Your Machine Configuration File

The machine configuration file is under `configs`. The configuration file contains the library paths, the compiler types, the compilation flags, and the GPU compute capability. Check the existing configuration files to see if there is one that matches your machine.

If the machine configuration file is not available for your machine or the existing one is not appropriate, you will need to create a new one. Follow the instructions in [[Machine Configuration File | Installation:-Machine-Configuration-File]] to set the configuration file of your machine.

If there is a configuration file that matches your machine, you should set it as the default by

```bash
sh tool/config/set_settings.sh --local --machine=your_machine
```

For example, setting `--machine=pleiades` with the above command will use the `configs/pleiades.config` machine configuration when compiling the code.

### Quick Demos

Read the following guides to learn how to configure, compile, and run GAMER.

1. [[1D Shock Tube: OpenMP with/without GPU acceleration | Quick-Start:-1D-Shock-Tube]]
2. [[3D Blast Wave: hybrid MPI/OpenMP/GPU + yt analysis | Quick-Start:-3D-Blast-Wave]]

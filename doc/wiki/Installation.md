1. Set up the machine configuration file

   Please see [[Machine Configuration File | Installation:-Machine-Configuration-File]].

2. Go to the source directory

   ```bash
   cd src
   ```

3. Generate `Makefile` using the Python script `configure.py`

   To get the `Makefile`, please execute the following command:

   ```bash
   python configure.py --machine=your_configuration_file [--your_arguments]
   ```

   `your_configuration_file` is the configuration filename you got from step 1, and `[--your_arguments]` should match your simulation requirements. Please check out [[Option List | Installation:-Option-List]] for all the available options.

   For example, the following command uses the `configs/pleiades.config` machine configuration, sets the FFTW method to `FFTW2`, and enables gravity and GPU.

   ``` bash
   python configure.py --machine=pleiades --fftw=FFTW2 --gravity=true --gpu=true
   ```

> [!TIP]
> An example script `generate_make.sh` to generate Makefile can be found in each test problem folder,
e.g., `example/test_problem/Hydro/AcousticWave/generate_make.sh`.

4. Compile the code

   ```bash
   make clean
   make
   ```

> [!TIP]
> To reduce the compilation time, you can perform a parallel
compilation by `make -j N`, where `N` is the number of compilation
jobs to run in parallel. For example, the following command will
invoke 4 compilation jobs simultaneously:
> ```bash
> make -j 4
> ```
> However, please consult the documentation of your system to avoid
violating the usage policy.

If the compilation succeeds, you will see the following message
<pre>
Compiling GAMER --> Successful!
</pre>
and get an executable `gamer`, which will be automatically copied to `../bin/gamer`.

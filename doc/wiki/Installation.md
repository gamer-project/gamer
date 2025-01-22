1. Set up the machine configuration file

   Please see [[Machine Configuration File | Installation:-Machine-Configuration-File]].

<a name="default_setting"></a>

2. Set your machine configuration file as default

   ```bash
   sh tool/config/set_settings.sh --local --machine=your_machine
   ```

> [!NOTE]
> If you want to set the default machine configuration file for all GAMER copies under your user account, use the `--global` option instead of `--local`.
You can still override the global setting for individual GAMER copies using the `--local` option.
Furthermore, you can override the default setting by passing the [[--machine | Installation:-Option-List#--machine]]=`your_machine` when executing `configure.py`. 

3. Go to the source directory

   ```bash
   cd src
   ```

4. Generate `Makefile` using the Python script `configure.py`

   To get the `Makefile`, please execute the following command:

   ```bash
   python configure.py [--your_arguments]
   ```

   `[--your_arguments]` represent the options that should align with your simulation requirements. Refer to [[Option List | Installation:-Option-List]] for a complete list of available options.

   For example, the following command sets the FFTW method to `FFTW2` and enables gravity and GPU.

   ``` bash
   python configure.py --fftw=FFTW2 --gravity=true --gpu=true
   ```

> [!TIP]
> An example script `generate_make.sh` to generate Makefile can be found in each test problem folder,
e.g., `example/test_problem/Hydro/AcousticWave/generate_make.sh`.

> [!TIP]
> Since there are too many options in GAMER, we introduce the autocomplete feature of `configure.py`. You can set the feature by the following steps:
> 1. Update the [shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)) of `configure.py` if needed
>
>    For example, replace `#!/usr/bin/python3` with your `python` path.
>    <details>
>    <summary><u><i>top of configure.py</i></u></summary>
>    <pre>
>    #!/usr/bin/python3
>    """
>    User guides of this script are provided in the following link.
>
>    https://github.com/gamer-project/gamer/wiki/Installation
>
>    """
>    </pre>
>    </details>
>
> 1. Copy the autocomplete shell script to your `/home/usr` (`~`)
>    ```bash
>    cp tool/bash/config_autocomplete.sh ~/
>    ```
>
> 1. Update the `~/.bashrc` to load the autocomplete script
>
>    Please add the following line to `~/.bashrc`:
>    ```bash
>    source ~/config_autocomplete.sh
>    ```
>    NOTE: There are three commands registered in the `config_autocomplete.sh`, `python` `python3` and `./configure.py`.
>          If you would like to register your command, please add the command at the bottom of the script.
>
> 1. Reload `~/.bashrc` to enable the feature
>    ```bash
>    source ~/.bashrc
>    ```
>
> Now, try to type `./configure.py` then press `<tab>` multiple times!

5. Compile the code

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

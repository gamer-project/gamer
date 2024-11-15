## Compile GAMER
1. Setup machine configuration file
   
   Please see [[Machine Configuration File | Installation:-Machine-Configuration-File]].

2. Go to the source directory

   ```bash
   cd src
   ```

3. Generate `Makefile` using the Python script `configure.py`

   Please check out [[Generate Makefile | Installation:-Generate-Makefile]].

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

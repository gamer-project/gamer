This test problem Plummer demonstrates the in situ Python analysis feature in GAMER.

***

1. Install the required packages for this demo.
  * [libyt and yt_libyt](https://yt-project.github.io/libyt/HowToInstall.html#how-to-install): the in situ analysis library and its yt frontend.
> [!IMPORTANT]
> This example assumes it is using libyt interactive mode. Please compile libyt in interactive mode.
  * [yt](https://yt-project.org/): the core analysis tool.
  * [[FFTW | Installation: External Libraries#fftw]]: needed in this test problem.
  * [[HDF5 | Installation: External Libraries#hdf5]]: this is optional since we are not going to output any data.

2. Generate `Makefile` with the following options:
   ```bash
   --gravity=true --particle=true --passive=2 --mpi=true --fftw=FFTW2 --libyt=true --libyt_interactive=true
   ```
   Also please set the paths `FFTW_PATH`, `MPI_PATH`, and `LIBYT_PATH` in [[configuration file | Installation:-Machine-Configuration-File ]].

3. Compile the code in the `src` folder.

   ```bash
   make clean
   make -j 4
   ```

<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
   ...
   ...
Compiling GAMER --> Successful!
</pre>
</details>

4. Create a working directory in the `bin` folder and copy the GAMER executable and Plummer test problem files.

   ```bash
   mkdir Plummer                                       # create a folder for this test
   cd Plummer
   cp -r ../../example/test_problem/Hydro/Plummer/* .  # copy test problem settings
   cp ../../src/gamer .                                # copy GAMER executable
   ```

5. In `Input__Parameter`, set `OPT__OUTPUT_TOTAL` to `0` since we don't need to dump data to disk in this demo.
(If you want to make GAMER dump data snapshots, enable `SUPPORT_HDF5` and set `HDF5_PATH` in `Makefile`. Keep `OPT__OUTPUT_TOTAL` to `1`.)

   ```
   # data dump
   OPT__OUTPUT_TOTAL       0     # output the simulation snapshot: (0=off, 1=HDF5, 2=C-binary) [1]
   ```

6. In `Input__TestProb`, set `Plummer_Collision` to `1` for simulating colliding Plummer clouds.

   ```
   Plummer_Collision       1     # (0/1)--> single Plummer cloud/two colliding Plummer clouds
   Plummer_AddColor        1     # assign different colors to different clouds (must turn on Plummer_Collision
                                 # and set NCOMP_PASSIVE_USER to 2 in the Makefile) [0]
   ```

7. Create a `LIBYT_STOP` file. This demo assumes libyt is in interactive mode.
It will only enter interactive Python prompt if errors occur while running Python codes
or it detects the file `LIBYT_STOP`.

   ```bash
   touch LIBYT_STOP
   ```

8. Run GAMER.

   ```bash
   OMPI_MCA_osc=sm,pt2pt mpirun -np 4 ./gamer
   ```

9. Results from the inline script.

Inline script `inline_script.py`:
```python
import yt_libyt
import yt

yt.enable_parallelism()

def yt_inline():
    # Get data
    ds = yt_libyt.libytDataset()

    # Do ProjectionPlot to field Cloud0.
    sz = yt.ProjectionPlot(ds, 'z', ('gamer', 'Cloud0'), center='c')

    # Do ParticlePlot
    par = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y', 'particle_mass', center='c')

    if yt.is_root():
        sz.save()
        par.save()

def yt_inline_inputArg( fields ):
    pass
```

Projection of field `Cloud01` along the z-axis:

[[/images/InSitu_Projection_z_Cloud0.png | alt=insitu_projection]]


Particle plot using particle position x and y as axes with particle mass mapped to a colorbar:

[[/images/InSitu_Particle_z_particle_mass.png | alt=insitu_particle]]

10. Since we have run libyt in interactive mode, the program should now pause and wait for user input.

<pre>
=====================================================================
  Inline Function                              Status         Run
---------------------------------------------------------------------
  * yt_inline                                  success         V
  * yt_inline_inputArg                         success         V
=====================================================================
>>>
</pre>

11. You can explore [interactive prompt](https://yt-project.github.io/libyt/InSituPythonAnalysis/InteractivePythonPrompt.html)
or simply press `Ctrl + C` to shut down the process.

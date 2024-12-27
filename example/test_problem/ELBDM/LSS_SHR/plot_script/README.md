# Analysis of the Soliton-Halo relation in the ELBDM model

This script is used to analyze the soliton-halo relation in the ELBDM model. The script reads the data from the simulation output stored at upper directories and plots the soliton-halo relation.

## Usage

1. Set the parameters in the script. The parameters are:
    - `file_in`: The path to the simulation output data.
    - `START_ID`: The starting ID of the simulation output.
    - `END_ID`: The ending ID of the simulation output.
    - `DELTA_ID`: The interval between two simulation outputs.
    - `HALO`: The halo ID.

    ```bash
    file_in=../Data
    START_ID=51
    END_ID=68
    DELTA_ID=1
    HALO=1
    ```

2. Create some directories to store the output data.
    ```bash
    mkdir -p prof_dens prof_mass volume velocity proj slice
    mkdir -p  velocity/v2_of_radius velocity/jeans "velocity/output_$HALO" 
    ```

3. Plot the projection of data in x and z directions to find the halo center.
    ```bash
    mpirun -np 4 python3 plot__proj-z-dens.py -s $START_ID -e $END_ID -halo $HALO -axis z
    mpirun -np 4 python3 plot__proj-z-dens.py -s $START_ID -e $END_ID -halo $HALO -axis x
    ```
    After running the above commands, estimate the halo center from the projection plots by eye.

4. Output the density profile and some properties of the halo using estimated halo center.
    Change the `coordinates = ts[0].arr( [ 0, 2.2, 0.4 ] , 'code_length' )`  in the script `plot__profile_mpi.py` to the estimated halo center.
    ```bash
    mpirun -np 4 python3 plot__profile_mpi.py -s $START_ID -e $END_ID -halo $HALO
    ```
    It will output the  profile at /prof_dens/, /prof_mass/, and /volume/ directories. The properties of the halo are stored in the file `Halo_Parameter_$HALO` in physics coordinate.

5. Plot the density profile of the halo.
    ```bash
    python3 plot__profile_with_data.py -s $START_ID -e $END_ID -halo $HALO
    ```
    The figures will be stored in the directory `prof_dens`.

6. Calculate the velocity profile of the halo.
    This would need `GAMER_EXTRACT_PROFILE` to calculate the velocity profile. Please see /tool/analysis/gamer_extract_profile and compile the code.
    The `MultiFiles__GAMER_ExtractProfile` can read `Halo_Parameter_$HALO` and run `GAMER_EXTRACT_PROFILE` to calculate the velocity profile of a sphere with 1.5 times radius halo.
    ```bash
    source MultiFiles__GAMER_ExtractProfile.sh $file_in $START_ID $END_ID $DELTA_ID $HALO
    mv Ave* Vir* velocity/output_$HALO
    ```
    The velocity profile will be stored in the directory `velocity`.

7. Plot the velocity profile of the halo.
    ```bash
    python3 plot__velocity.py -s $START_ID -e $END_ID -halo $HALO
    ```
    This would plot the velocity profile and calculate the Jeans velocity. (Jeans velocity would store in the directory `velocity/jeans` and would not calculate if the files already exists.) The figures will be stored in the directory `velocity`. It also calculate the average velocity and  output the result in the file `halo_velocity_$HALO`.

8. Plot the ratio of factors.
    ```bash
    python3 plot__all.py -s $START_ID -e $END_ID -halo $HALO
    ```

9. You can also plot the slice of the halo. It will plot the grid so that you can see if the resolution is enough.
    ```bash
    mpirun -np 4 python3 plot__slice.py -s $START_ID -e $END_ID -halo $HALO
    ```
    The figures will be stored in the directory `slice`.

10. The script that reproduce the figure in the paper is in the directory `/script_for_paper_fig`. You need to finish the above steps to run the script.
    - `load_simulation.py` is used to load the simulation data. Use `load(path, halo, start_id, end_id, name)` to load the data.
    - `SHR_compare_paper.py` is used to reproduce the figure. 1 in the paper.
    - `plot__profile_paper.py` is used to reproduce the figure. 2 in the paper.
    - `plot__velocity_paper.py` is used to reproduce the figure. 3 in the paper.
    - `plot__c_nonisothermality_paper.py` is used to reproduce the figure. 4 in the paper.
    - `plot__ratio_paper.py` is used to reproduce the figure. 5 in the paper.

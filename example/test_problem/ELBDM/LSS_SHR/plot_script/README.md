# Analysis of the Soliton-Halo relation in the ELBDM model

These scripts are used to analyze the soliton-halo relation in the ELBDM model. The analysis includes generating projections, density profiles, velocity profiles, and other properties from simulation outputs.

## Features

- Projection plotting to estimate halo centers.
- Density and velocity profile extraction.
- Calculation of Jeans velocity and other physical properties.
- Scripts to reproduce figures from related research papers.

## Requirements

1. **Software Requirements**:
- Python 3
- MPI for parallel execution
- GAMER simulation tools (`GAMER_EXTRACT_PROFILE`)

2. **Dependencies**:
Ensure required Python libraries are installed, such as `numpy`, `matplotlib`, `pandas`, `scipy`, `yt`, etc.

3. **Simulation Data**:Place your simulation output data in a directory, e.g., `../Data`.


## Usage

### 1. Set Parameters

Edit the script parameters before execution:

```bash
file_in=../Data # Path to the simulation output data
START_ID=51 # Starting ID of the simulation output
END_ID=68 # Ending ID of the simulation output
DELTA_ID=1 # Interval between two simulation outputs
HALO=1 # Halo ID
```

### 2. Create Output Directories

```bash
mkdir -p prof_dens prof_mass volume velocity proj slice
mkdir -p  velocity/v2_of_radius velocity/jeans "velocity/output_$HALO" 
```

### 3. Plot Data Projections

Generate x and z projections to estimate the halo center:

```bash
mpirun -np 4 python3 plot__proj-z-dens.py -s $START_ID -e $END_ID -halo $HALO -axis z
mpirun -np 4 python3 plot__proj-z-dens.py -s $START_ID -e $END_ID -halo $HALO -axis x
```
Outputs are saved in `proj/`. Manually estimate the halo center from the projection plots.

### 4. Extract Density Profile

Update the coordinates variable in plot__profile_mpi.py with the estimated halo center:
```python
coordinates = ts[0].arr([0, 2.2, 0.4], 'code_length')
```
Run the script:
```bash
mpirun -np 4 python3 plot__profile_mpi.py -s $START_ID -e $END_ID -halo $HALO
```
Outputs:
- Density profile in `prof_dens/`.
- Halo properties in `Halo_Parameter_$HALO`.

After calculating the halo properties, the script automatically sets the calculated center of the current halo as the estimated center for the next halo.

### 5.  Plot Density Profile

Generate density profile plots:
```bash
python3 plot__profile_with_data.py -s $START_ID -e $END_ID -halo $HALO
```
Figures are saved in `prof_dens/`.

### 6. Calculate Velocity Profile

Compile GAMER_EXTRACT_PROFILE from `/tool/analysis/gamer_extract_profile` and use it to calculate the velocity profile:
```bash
source MultiFiles__GAMER_ExtractProfile.sh $file_in $START_ID $END_ID $DELTA_ID $HALO
mv Ave* Vir* velocity/output_$HALO
```
Outputs are saved in `velocity/`.

### 7. Plot Velocity Profile
Generate velocity profile plots and calculate Jeans velocity:
```bash
python3 plot__velocity.py -s $START_ID -e $END_ID -halo $HALO
```
Outputs:
- Velocity profiles in `velocity/`.
- Jeans velocity in `velocity/jeans`.
- Average velocity in `halo_velocity_$HALO`.

### 8. Analyze and Plot Ratios
Plot ratio factors:
```bash
python3 plot__all.py -s $START_ID -e $END_ID -halo $HALO
```

###9. Visualize Halo Slices
Generate slice plots to check grid resolution:
```bash
mpirun -np 4 python3 plot__slice.py -s $START_ID -e $END_ID -halo $HALO
```
Figures are saved in slice/.

### 10. Reproduce Research Figures
Use the scripts in /script_for_paper_fig to reproduce figures from the associated paper:
- `load_simulation.py`: Load simulation data. Change the file path and halo ID to load the simulation data.
  ```python
  load(path, halo, start_id, end_id, name)
  ```
- `SHR_compare_paper.py`: Reproduce Figure 1.
- `plot__profile_paper.py`: Reproduce Figure 2.
- `plot__velocity_paper.py`: Reproduce Figure 3.
- `plot__c_nonisothermality_paper.py`: Reproduce Figure 4.
- `plot__ratio_paper.py`: Reproduce Figure 5.

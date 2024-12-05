# Spectral Interpolation in GAMER
This section provides an overview of the spectral interpolation method in GAMER implemented via the Gram-Fourier continuation algorithm in the Jupyter notebook `GramFE.ipynb` located in `tool/table_maker/GramFE`. This algorithm is essential for computing extension and interpolation tables.

## Setup and Configuration

### Compilation and Runtime Settings
- **Compile-Time Flag**: Ensure GAMER is compiled with `SUPPORT_SPECTRAL_INT` (or configured with `--spectral_interpolation`).
- **Runtime Parameters**: Set `OPT__FLU_INT_SCHEME` and `OPT__REF_FLU_INT_SCHEME` to `8` for enabling spectral interpolation. Set `SPEC_INT_TABLE_PATH` to the directory containing `interpolation_tables` and `boundary2extension_tables`. Enable `SPEC_INT_XY_INSTEAD_DEPHA` to interpolate real and imaginary parts instead of density and phase around vortices, which has the advantage of being well-defined across vortices. `SPEC_INT_VORTEX_THRESHOLD` sets the vortex detection threshold for `SPEC_INT_XY_INSTEAD_DEPHA`. If the laplacian of the phase field exceeds this threshold, we interpolate the real and imaginary parts. `SPEC_INT_GHOST_BOUNDARY` sets the ghost boundary size for spectral interpolation. A large ghost boundary increases interpolation accuracy but can negatively affect performance.

### Obtaining Interpolation Tables

#### Downloading Interpolation Tables

- **Script**: Use `download_spectral_interpolation_tables.sh` located in `example/test_problem/ELBDM/LSS_Hybrid`.

#### Generating Interpolation Tables

- **Script**: Use `compute_interpolation_tables.py` located in `tool/table_maker/GramFE`.
- **Execution**: Run `mpirun -n 16 python3 compute_interpolation_tables.py` for table generation.
- **Time Consumption**: The script may take several hours to execute without additional output.

## Advantages and Disadvantages of the Spectral Interpolation Algorithm

### Advantages
- **Spectral Accuracy**: Offers high-precision interpolation results.
- **Flexibility**: Suitable for various interpolation needs within the GAMER framework.

### Disadvantages
- **Non-conservative and Non-monotonic**: Interpolation results may not always preserve these properties.
- **Computationally Intensive**: Spectral interpolation is slower than polynomial interpolation with a local stencil (i.e. interpolation options `1=MinMod-3D, 2=MinMod-1D, 3=vanLeer, 4=CQuad, 5=Quad, 6=CQuar, 7=Quar`).

## Testing the Spectral Interpolation Algorithm
Different algorithms in GAMER introduce different discretisation errors:
- numerical PDE solvers
- interpolation (new patches during refinement and ghost boundaries for neighbouring patches on different AMR levels)
- time interpolation for flexible time-steps (first-order)
- restriction operation (first-order)
- flux-fixup operation to ensure mass conservation

### Minimising errors

To minimise these errors and measure the performance of different interpolation algorithms, you may consider the following settings
- Use vortex pair test problems
- **Compile time parameters**: `ELBDM_SCHEME=ELBDM_WAVE`, `WAVE_SCHEME=WAVE_GRAMFE` and `GRAMFE_SCHEME == GRAMFE_MAMTUL`
- **Runtime parameters**: Turn `OPT_FIXUP_FLUX`, `OPT_INIT_RESTRICT`, `OPT_FIXUP_RESTRICT`, `OPT_PHASE_INT` off, set `OPT__DT_LEVEL = 1`

## Known Issues

### Negative Density in Spectral Interpolation

- **Problem Description**: The spectral interpolation method may introduce negative density values. This issue arises because the interpolation is non-monotonic.
- **Consequences**: Negative density can lead to NaN (Not a Number) values in both density and flux calculations within the fluid scheme. This problem has been observed to potentially cause `MaxPot == 0.0` errors in certain simulations.
- **Current Investigation**: The issue is under active investigation, with more details and ongoing discussions available on GitHub and Slack. Specific instances of this problem have been shared [here](https://calab-ntu.slack.com/archives/C02JV33GF5Y/p1693032184198239).
- **Proposed Solution**: A preliminary solution involves applying a density floor when preparing density ghost zones in `Flu_Prepare()`. This adjustment could prevent the generation of negative density values during spectral interpolation. Further details and implementation suggestions can be found in the related GitHub comment thread.

#### Recommendations for Users

- Users encountering similar issues or observing unexpected results in density and flux calculations should use GAMER's other more reliable interpolation algorithms (i.e. `6=CQuar` (conservative quartic interpolation)).


## Implementation Details of the Spectral Interpolation Algorithm

### Overview
- The spectral interpolation method in GAMER is implemented through a Gram-Fourier continuation algorithm in the Jupyter notebook `GramFE.ipynb`, located in `tool/table_maker/GramFE`. The basic idea is to find a periodic extension of a given interpolant. This periodic function can then be interpolated using a DFT. This interpolation is not monotonic and not conservative, but highly accurate.
### Algorithm References
- Lyon and Bruno's Original Work:
  - [Journal of Computational Physics, 2009](https://doi.org/10.1016/j.jcp.2009.11.020)
  - [Journal of Computational Physics, 2010](https://doi.org/10.1016/j.jcp.2010.01.006)

### Key Components and Libraries
- **Arbitrary Precision Arithmetic**: Utilizes `mpmath` library.
- **Gram-Schmidt Orthogonalization**: Implements this algorithm for basis computation.
- **SVD Algorithm**: Used for high-precision continuations of Gram polynomials.
- **FFT Computations**: Involves forward and backward FFTs in high precision.
### Algorithm

The Python code generates two types of tables essential for the spectral interpolation algorithm, stored in `interpolation_tables` and `boundary2extension_tables`.

#### Interpolation for Small Input Sizes
- **Precomputed Operation**: The entire interpolation process (`GramFE extension * DFT * shift in k space * IDFT`) is precomputed as a single matrix.
- **Storage**: The resulting matrices are stored in double-precision binary files named `N.bin`, where `N` is the input interpolation size.
- **Efficiency**: This approach is efficient for small values of `N`, requiring `N^2` operations for each matrix multiplication at runtime.
- **Equivalent to Polynomial Interpolation**: Notably, for `m = N`, this method is equivalent to polynomial interpolation of order `N`, as confirmed by numerical tests.

#### Interpolation for Large Input Sizes
- **Runtime Computation**: For large inputs, the extension and DFTs are computed at runtime using the FFT algorithm.
- **Table Utilization**: The tables `%ld_%ld_%d_%d.bin` are used, where `nd`, `m`, `Gamma`, `g` represent the table parameters. By applying these tables to the boundary data (size `2 * m`), an extension of size `nd` is produced.
- **Adaptive Size Selection**: The size `nd` is adaptively chosen during runtime to optimize the FFT algorithm's performance, favoring sizes with small prime factorizations.
- **Matrix and FFT Operations**: The interpolation involves `2 * m * m` matrix operations and additional `O(N log N)` operations for the FFT, managed by the FFTW library.

#### Accuracy Considerations
- **Determining Factors**: The interpolation's accuracy heavily depends on the quality of the periodic extension.
- **Optimal Parameters**: Parameters `Gamma = 150`, `g = 63`, `nd = 32`, and `m = 8` strike a balance between stability and accuracy.
- **Artifacts and Behavior**: Higher values of `m` can lead to interpolation artifacts, while lower values might reduce accuracy but result in smoother interpolation.

The following plots show a number of accuracy comparisons between quartic interpolation and spectral interpolation:

![sine](https://github.com/hyschive/gamer-fork/assets/6187378/34e2f303-916b-443a-8e18-097adfbe32ff)
![exp](https://github.com/hyschive/gamer-fork/assets/6187378/09ec67cb-3362-479c-a901-d212011fcc24)
![high-freq](https://github.com/hyschive/gamer-fork/assets/6187378/374cf84f-4dc1-4a09-86f8-47b31c23b05c)
![runge-function](https://github.com/gamer-project/gamer/assets/6187378/837c1519-1a5b-4425-95f5-8d8f4f9f599b)

The last function demonstrates that while the high-order boundary polynomials used in the boundary matching of the GramFE algorithm suffer from the [Runge phenomenon](https://en.wikipedia.org/wiki/Runge%27s_phenomenon), for larger N the numerical properties of the algorithm are determined by the plane wave basis and convergence sets in.
Potentially, the large error of the GramFE method in this case can be remedied by switching to lower-order boundary polynomials for small N at the cost of sacrificing accuracy for more well-behaved functions.



## Performance
The following plots show a performance comparison between quartic interpolation and spectral interpolation using matrix multiplication and FFT. It demonstrates that for $N < 30$, the matrix multiplication interpolation is faster.

![timing](https://github.com/hyschive/gamer-fork/assets/6187378/1e2f5ecc-8003-45ac-9886-47bb1811ea31)


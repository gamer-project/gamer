# Spectral Solver for the Schrödinger Equation in GAMER

## Setup and Configuration

### Compilation Options
- **Wave Scheme Selection**: To enable the spectral solver, set `WAVE_SCHEME` to `WAVE_GRAMFE` in the Makefile or configure with `--wave_scheme=GRAMFE` for the ELBDM model.
- **Spectral Scheme Options**: Choose between `MATMUL` (faster for `PATCH_SIZE=8`) and `FFT` (faster for larger patch sizes) by setting `SCHEME` in the Makefile or configure with `--gramfe_scheme` accordingly.
- **Additional Libraries**: The CPU version requires FFTW 3 for single precision and FFTW 2/3 for double precision. The GPU version of `FFT` needs the cuFFTDx library. Set `FFTW2_PATH`, `FFTW3_PATH`, and `CUFFTDX_PATH` in the Makefile or configuration files `*.config`.

### Precision and Performance Options
- **Floating Point Precision**: Controlled by `GRAMFE_ENABLE_SINGLE_PRECISION`. The default is double precision. Single precision is not recommended due to stability issues.
- **Extension Parameters**: Set via `GRAMFE_ORDER`, `GRAMFE_NDELTA`, and `GRAMFE_ND` in `Macro.h`. New extensions can be added in `GramExtensionTables.h`.

## Advantages and Disadvantages

### Advantages
- **Higher Accuracy**: Achieves greater accuracy than the inbuilt finite-difference scheme.
- **Efficient FFT Performance**: Optimized for factorizations into small prime numbers.

### Disadvantages
- **Computational Cost**: Slightly larger computational cost compared to the finite-difference scheme.
- **GPU Solver Limitations**: GPU solver is slower on non-HPC GPUs for required double-precision operations.
- **No Mass Conservation**: Scheme is inherently not conservative. (Mass can move from physical into extended domain.)
- **Stability only Empirical**: Requires filter because of spectral blocking.

## Test problems
- Observe differences in accuracy for all inbuilt test problems that output L1 error
- Check test problem `RestrictionMismatch` to see how spectral solver can reduce artifacts
     - It features a traveling Gaussian wave refined around the center.
     - The following plots generated with the default visualisation scripts of the test problem highlight the mismatch between the wave function on level 0 and level 1 before and after time evolution.
     - Using the finite-difference scheme on all levels, there is a large mismatch between the evolution of the wave function on different levels.
![FD1_1](https://github.com/gamer-project/gamer/assets/6187378/22fdc630-3caf-4134-beec-0fa9f8efcc78)
![FD1_2](https://github.com/gamer-project/gamer/assets/6187378/3385ddaf-0c19-4525-9674-c27049958c1a)
     - Using the base-level spectral method together with the finite-difference scheme on level 1 reduces the mismatch slightly.
![FD2_1](https://github.com/gamer-project/gamer/assets/6187378/40652423-025f-4376-b227-736134a71dff)
![FD2_2](https://github.com/gamer-project/gamer/assets/6187378/c55eda48-ddc0-46b5-847a-006099d4033a)
     - Using the base-level spectral method together with the local spectral scheme on level 1 significantly reduces artifacts.
![FD3_1](https://github.com/gamer-project/gamer/assets/6187378/833975c8-54eb-40b7-9d4b-7cddbf994ed9)
![FD3_2](https://github.com/gamer-project/gamer/assets/6187378/65fc7c65-84ed-4859-baf0-9c3ad95f8d22)

## Implementation Details

### Overview
- The spectral solver is based on the Gram-Fourier extension algorithm, offering a local spectral approach to the Schrödinger equation.

### Gram-Fourier extension algorithm
- Given a non-periodic function on $[0, L]$, how to find a periodic extension on $[0, L + d]$?
- Optimise linear, under-determined least squares problem via SVD $A = U \Sigma V^T$.
![OptimisationProblem](https://github.com/gamer-project/gamer/assets/6187378/bf8ecfc5-174c-43f6-8684-2fbcfb462dc1)
- This approach gives excellent periodic extension up to arbitrary precision with the right choice of hyperparameters.
![FourierExtension](https://github.com/gamer-project/gamer/assets/6187378/4c55a306-a573-4966-aa11-14311f64f20e)
(Figure taken from [Mark Lyon's PhD thesis](https://thesis.library.caltech.edu/2992/), DOI: 10.7907/3FFW-GK56)
- **Drawback**: Cost for SVD of $n\times m$ matrix $\mathcal{O}(\min(m, n)^2 \cdot \max(m, n))$.
- How to make SVD extensions efficient?
- Key idea: Use **linearity** of optimisation problem.
- Expand a non-periodic function in terms of a suitable basis at the left and right boundaries.
- Compute SVD extensions (orange) of basis functions (blue) (left: constant polynomials, center: linear polynomial, right: quadratic polynomial).
![Extensions](https://github.com/gamer-project/gamer/assets/6187378/2662f3ef-dd0b-4d84-8a06-d8cf41d4b2a0)
- Express periodic extension as a suitable linear combination of SVD basis functions at the left and right boundaries.
![Boundary](https://github.com/gamer-project/gamer/assets/6187378/2f4faa42-0aa0-4057-aba8-67d0ef1a5e24)
![ApplyExtension](https://github.com/gamer-project/gamer/assets/6187378/5c0acbb1-b77d-40ba-ad07-f8b9f052cf12)
(Figures taken from [Mark Lyon's PhD thesis](https://thesis.library.caltech.edu/2992/), DOI: 10.7907/3FFW-GK56)

### Steps of the Algorithm
1. **Boundary Polynomial Expansion**: Expand real and imaginary parts at boundaries (not the same as GAMER's ghost boundaries).
2. **Periodic Extension**: Continue polynomials using precomputed periodic extensions. This extension is unphysical and doesn't use neighboring patch information.
3. **FFT Evolution**: Use FFT on the extended domain to evolve the wave function, which becomes periodic on this domain.
4. **Avoidance of Gibbs Phenomenon**: Smooth extensions largely prevent the Gibbs phenomenon.
5. **Discarding Extended Domain**: After computation, the extended domain is discarded.

## Differences between MATMUL and FFT Schemes

### Overview
Understanding the differences between the `MATMUL` and `FFT` schemes in GAMER is crucial for users to choose the appropriate method for their specific simulation needs. Both schemes are based on the Gram-Fourier extension algorithm but differ in their computational approach and performance characteristics.

### MATMUL Scheme

#### Description
- **Matrix Multiplication Approach**: In the `MATMUL` scheme, the entire operation sequence - GramFE extension, FFT, time evolution, IFFT, and discarding of extension data - is precomputed and represented as a single matrix (which is possible because of linearity).
- **Application to Input Vectors**: This precomputed matrix can be directly applied to input vectors, eliminating the need to compute each operation at runtime.

#### Advantages
- **Reduced Runtime Computation**: By precomputing the entire operation sequence, the `MATMUL` scheme minimizes runtime computational overhead.
- **Efficiency for Small Patch Sizes**: This scheme is particularly efficient for smaller data sizes where the cost of matrix multiplication is lower.

#### Disadvantages
- **Memory Usage**: The precomputed matrix, especially for large patch sizes, can consume significant memory resources.
- **Precomputation Time**: There is an upfront cost in time to compute and store the matrix, which might be substantial depending on the accuracy used.

### FFT Scheme

#### Description
- **Runtime Computation**: The `FFT` scheme performs the GramFE extension, FFT, time evolution, IFFT, and discarding of extension data operations at runtime.

#### Advantages
- **Scalability for Large Patch Sizes**: This scheme is more scalable for larger patch sizes, as it does not require storing large precomputed matrices.
- **Better Round-off Errors**: Slightly better mass conservation properties than matrix multiplication scheme. This gap may be overcome by using higher precision for the precomputation of the time evolution matrix in the matrix multiplication scheme.
#### Disadvantages
- **Increased Runtime Operations**: Each step of the operation sequence is computed at runtime, which can increase the computational overhead, especially for larger simulations.
- **Dependence on FFT Performance**: The efficiency of the scheme is closely tied to the performance of the FFT algorithm used.

### Conclusion
- **Choice of Scheme**: By default `MATMUL` should be used. `FFT` should be considered for larger patch sizes and debugging.

## Performance Metrics
- The performance ratio `Perf_FD / Perf_Spectral` indicates the efficiency compared to the finite-difference scheme. This varies based on patch size (e.g., 4.3 for `PATCH_SIZE=8` and 2.6 for `PATCH_SIZE=16` for `FFT`; comparable performance of `MATMUL` and `FD` for `PATCH_SIZE=8`).

## Issues and Considerations
- **Stability with Single Precision**: Using single precision can lead to instability in the scheme.
- **Patch Size Impact**: The choice of patch size directly affects the performance and efficiency of the FFT operations.
- **Not conservative**: Scheme is not conservative and there are no fluxes for the flux-fixup operation.

## Possible improvements
- **Stability and Accuracy**: Fine-tune the parameters of the GramFE extension, the boundary size as well as the filter parameters to achieve a good compromise between stability, accuracy and the size of the ghost zone.
- **Mass Conservation**: Solve continuity equation (possibly using Gram-FE and DFT) to ensure mass conservation and compute estimate for fluxes.
- **Density & Phase**: Compute extension for density & phase or another set of more suitable variables (e.g. density * sin(phase/N) and density * cos(phase/N) with N>1 to avoid phase discontinuity at vortices) for higher accuracy.

### High Precision Requirements for Evolution Matrix Calculation

#### Current Precision and Limitations

The calculation of the evolution matrix in the `MATMUL` scheme currently requires at least 128-bit precision to maintain an error below single precision. The extension matrix has a high condition number, typically around 10^15 or higher, which necessitates this high precision.

#### Proposed Enhancements

- **Upgrade to 256-bit Type**: Ideally, the `gramfe_evo_float` data type should be upgraded to a 256-bit type in future iterations. This would significantly improve the accuracy of the evolution matrix calculation.
- **Impact on Extension Tables**: Correspondingly, the `GramFE_ExtensionTables` would need to be updated to include more significant digits to leverage the increased precision effectively.
- **Current Usage of Long Double Type**: Some routines currently use the `long double` type instead of the `quadmath` version. This choice is based on the requirement that these routines do not need more than double precision. In light of the proposed shift to a 256-bit type for `gramfe_evo_float`, a review of the data types used in these routines may be warranted to ensure consistency and optimal performance across the scheme.

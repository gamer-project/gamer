* #### Running out of GPU memory
   * **Description**:
   ```
   ERROR : CUDA ERROR : out of memory !!
   ```
   * **Solution**: Reduce [[NSTREAM | GPU#NSTREAM]] until the issue is resolved.
* * *

* #### Instability with MHD
   * **Description**: See this [issue](https://github.com/gamer-project/gamer/issues/80#issuecomment-1030177067) report.
   * **Solution**: Update `CUDA` to 11.3 or higher.

* * *
* #### Check static arrays with [AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer) (ASan)
   * **Description**: Detect incorrect usage of static arrays.
   * **Steps**:
      1. Use the `g++` compiler.
      2. Compile and link with flags `-fsanitize=undefined -fsanitize=address`.

   * **Example (on the `eureka` machine)**
      1. Uncomment the following lines in `configs/eureka_gnu.config`
         ```
         #CXXFLAG -fsanitize=undefined -fsanitize=address
         #LIBFLAG -fsanitize=undefined -fsanitize=address
         ```
      2. Generate `Makefile`
         ```
         python configure.py --machine=eureka_gnu [--your_other_arguments]
         ```
      3. Compile
         ```
         make clean && make -j4
         ```
* * *

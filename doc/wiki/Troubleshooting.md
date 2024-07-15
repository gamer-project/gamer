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
* #### Check static array by [AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer) (ASen)
   * **Description**: To detect the wrong usage of the static array.
   * **Steps**:
      1. Use `g++` compiler.
      1. Compile and link with flags `-fsanitize=undefined -fsanitize=address`.

   * **Example of using `eureka` machine**
      1. Uncomment the following lines in `configs/eureka_gnu.config`
         ```
         #CXXFLAG -fsanitize=undefined -fsanitize=address
         #LIBFLAG -fsanitize=undefined -fsanitize=address
         ```
      1. Generate `Makefile`
         ```
         python configure.py --machine=eureka_gnu [--your_other_arguments]
         ```
      1. Compile
         ```
         make clean && make -j4
         ```
* * *

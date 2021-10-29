

****************************************************************************
GAMER_Uniform2D2gnuplot : convert 2D uniform data to the gnuplot pm3d format
****************************************************************************


Supported 2D input data :
============================================================================
1. Outputs of "GAMER" with "OPT__OUTPUT_PART=1/2/3" and "OPT__OUTPUT_BASE=1"
2. Outputs of "GAMER_GetCube" with the option "-n 1/2/3/4/5/6"

***In other words, it only works with "2D uniform data"***



Procedure :
============================================================================
1. Compile GAMER_Uniform2D2gnuplot.cpp by "make"

2. sh CreateFigures.sh InputName [StartID EndID]

3. A very simple script "CreateMovie.sh" can be used to create a GIF animation
   named "Movie.gif".
   --> Please refer to the ImageMagick website for more details.



Note :
============================================================================
1. Note that if one wants to transform multiple files automatically by
   specifying "StartID and EndID", the data ID suffix _XXXXXX should NOT be
   included in InputName.

   For example, to create figures for the three input files named
   "Input_000001, Input_000002, Input_000003", one should use

      sh CreateFigures.sh Input 1 3

   Three PNG images files "Input_000001.png, Input_000002.png,
   Input_000003.png" will be created.

2. If neither StartID nor EndID is provided, the InputName should contain
   the full file name.

3. The InputName can contain the absolute path of the targeted file. The
   absolute path will NOT be included in the output file name.

4. One does not need to provide the number of data columns and the rows.
   They will be determined from the input file automatically.

5. In the script "CreateFigures.sh", one can set the
   variable "GPT_TERM" to different gnuplot terminals to generate different
   image formats (e.g., postscript/emf/gif ... ). Type "help terminal" in
   gnuplot to see all terminals supported in gnuplot.

6. To modify any gnuplot setting, please modify the small gnuplot script
   "gnuplot_pm3d.gpt"

   Examples:
   (a) Plot pressure   --> splot "`echo $FILE_TEMP`" u 1:2:8
       (3/4/5/6/7/8/9 --> Density/Momentum x/y/z/Energy/Pressure/Potential)
   (b) No contour plot --> mark the line "set contour base"
   (c) Change color    --> set palette rgbformulae ?,?,?
   (d) Set x/y range   --> set xrange [?:?]
                           set yrange [?:?]

7. Use the option "-n" to set the number of header lines (default = 1)



Revision History :
============================================================================

Version 1.1.0  09/21/2011
-------------------------
1. Support loading arbitrary number of data columns
   --> The option "-p" is removed
2. To load a single file, one can provide the absolute file name without
   specifying the range of data IDs
3. More error checks in "CreateFigures.sh"
4. One can specify the number of header lines through the option "-n"


Version 1.0.1  02/11/2011
-------------------------
1. Support gnuplot version 4.0 or higher
   --> no longer pass the file information as gnuplot command-line arguments
   --> use environment variables (and retrieve them in gnuplot by `echo $VAR`)


Version 1.0    02/09/2011
-------------------------
1. First version


Bug :          02/09/2011
-------------------------


Unfinished :   02/09/2011
-------------------------





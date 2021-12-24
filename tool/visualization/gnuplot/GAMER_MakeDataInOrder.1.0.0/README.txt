

****************************************************************************
GAMER_MakeDataInOrder : Sort the input data into ascending numerical order
****************************************************************************



Procedure :
============================================================================
1. Compile program by "make"

2. sh MultiFiles__GAMER_MakeDataInOrder.sh InputName [StartID EndID [DeltaID]]



Note :
============================================================================
1. Note that if one can use the script MultiFiles__XXX.sh to work on multiple
   files automatically. Please specify "InputName StartID EndID [DeltaID]".
   Note that the data ID suffix _XXXXXX should NOT be included in InputName
   if StartID and EndID are provided.

   For example, to work on three input files named
   "Input_000001, Input_000002, Input_000003", one should use

      sh script MultiFiles__XXX.sh Input 1 3

2. If neither StartID nor EndID is provided, the InputName should contain
   the full file name.

3. If no output file name is specified, it will be set to
   "$InputName_Ordered" by default.

4. One does not need to provide the number of data columns and the rows.
   They will be determined from the input file automatically.

5. Use the option "-n" to set the number of header lines [default = 1].

6. Use the option "-c" to specify the targeded data column (the key) to be
   sorted (start with zero) [default = 0]



Revision History :
============================================================================

Version 1.0.0  11/29/2011
-------------------------
1. First version



Bug :          11/29/2011
-------------------------



Unfinished :   11/29/2011
-------------------------





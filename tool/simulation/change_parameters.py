#!/bin/python3
"""
This is a script of changing the parameters of Input__* files.

How to use it:
  1. Set up the files and the parameters you want to change.
    a. The file is set to be a `File()` class. The `File` takes four inputs: file name, constant
       parameters, changing parameters, and flag file or not. Please check out the `Classes` section
       for the detail of inputs.
    b. Wrap all the `File` classes to a list called `files`.

    * NOTICE:
      Parameter files:
        The first column specifies the option's name and the second column should be the value of the option.
      Flag files:
        The first line of a flag file should always start with `# ` and the following column
        names in the header specify the options' names.

  2. Overwrite the `execution()` function to your execution.
     For instance, `subprocess.called([kwargs["exe"]]).

  3. [Optional] To pass extra arguments to `execution()` from the command line,
     add your own `parser.add_argument` in the `Main` section.

An example of the script (from Alexander Kunkel):
-----------------------------------------------------------------------------------------------------
# Additional imports for interacting with file system
import shutil
import glob

...

def execution( **kwargs ):
    cwd               = os.getcwd()
    record_parameters = kwargs["record"]
    exe               = kwargs["exe"]


    # Run GAMER with executable passed as command line parameter
    subprocess.run([kwargs["exe"]], capture_output=False)

    # Analyse run: Call python scripts etc.
    os.system("python3 plot_wave_slice.py -s 1 -e 1")
    os.system("python3 analyse_l1_error.py")
    os.system("ls -alt > Record__FileSizes")

    # Create folder recording which parameters where changed
    par_dir  = exe[2:] + "_" + "_".join(record_parameters.keys())

    try:
        os.mkdir(par_dir)
    except FileExistsError:
        pass

    # Create folder recording current runtime parameters
    sub_dir  = str(record_parameters)
    dest_dir = os.path.join(par_dir, sub_dir)
    try:
        os.mkdir(dest_dir)
    except FileExistsError:
        pass

    # Move and or output files to folder, delete with os.remove() (files) and shutil.rmtree() (directories) if necessary
    for file in glob.glob(r'*.png'):
        shutil.move(os.path.join(cwd, file), os.path.join(cwd, dest_dir))

    for file in glob.glob(r'Input*'):
        shutil.copy(os.path.join(cwd, file), os.path.join(cwd, dest_dir))

    for file in glob.glob(r'Record*'):
        shutil.move(os.path.join(cwd, file), os.path.join(cwd, dest_dir))

    for file in glob.glob(r'VortexPairLinear*'):
        shutil.move(os.path.join(cwd, file), os.path.join(cwd, dest_dir))

    for file in glob.glob(r'Data*'):
        shutil.move(os.path.join(cwd, file), os.path.join(cwd, dest_dir))

    return

...

file_name1       = "Input__Parameter"
loop_paras1      = { "OPT__FIXUP_RESTRICT": [0, 1], "OPT__REF_FLU_INT_SCHEME": [4, 5, 6, 7],
const_paras1     = {"MAX_LEVEL": 1, "NX0_TOT_X": 64, "NX0_TOT_Y": 64, "OPT__INIT": 2, "END_STEP": 1000, "REFINE_NLEVEL": 1, "REGRID_COUNT":1, "OPT__INT_PHASE": 0}

file_name2       = "Input__Flag_Rho"
loop_paras2      = { "Density": [0, 1e3]}
const_paras2     = {}


file1  = File( file_name1, const_paras1, loop_paras1, False  )
file2  = File( file_name2, const_paras2, loop_paras2, True  )

files  = [ file1, file2 ]
-----------------------------------------------------------------------------------------------------

For developer:
1. The main concept is to use a recursion function instead of nest for loops, so the code stays clean
   and easy to maintain.
2. First we iterate the files then the parameters of each files. At the end of iteration, we called
   `execution()`.
"""
#====================================================================================================
# Import packages
#====================================================================================================
import numpy as np
import argparse
import re
import os
import sys
import subprocess



#====================================================================================================
# Global variables
#====================================================================================================
RECURSION_LIMIT = 1000  # the recursion depth limit (default in Python is 1000)
RETURN_FAIL     = False
RETURN_SUCCESS  = True



#====================================================================================================
# Classes
#====================================================================================================
class File():
    def __init__( self, file_name, consts, paras, flag_file ):
        """
        file_name : string. The file name.
        consts    : dict. The const parameter to be specified.
        paras     : dict with list element. The parameter to be chenged.
        flag_file : bool. The flag file or not.
        """
        self.name   = file_name
        self.paras  = paras
        self.consts = consts
        self.flag   = flag_file

        self.set_constants()

    def set_constants( self ):
        replace_func = replace_parameter_flag if self.flag else replace_parameter
        for key, val in self.consts.items():
            replace_func( self.name, key, val )
        return



#====================================================================================================
# Functions
#====================================================================================================
def iter_files( files, record_changed={}, **kwargs ):
    if len(files) == 0: return execution( record=record_changed, **kwargs )

    files_copy = list(files)
    f_class    = files_copy.pop(0)
    return iter_file_parameters( f_class.name, f_class.paras, f_class.flag, record_changed, files_copy, **kwargs )

def iter_file_parameters( file_name, paras, flag_file, record, rest_files=[], **kwargs ):
    if len(paras) == 0: return iter_files( rest_files, **kwargs )

    kwargs_copy  = kwargs.copy()
    paras_copy   = paras.copy()
    replace_key  = next( iter(paras_copy) ) # take the first key to replace
    replace_func = replace_parameter_flag if flag_file else replace_parameter
    vals         = paras_copy.pop(replace_key)

    for val in vals:
        if not kwargs["quite"]: print("File %-25s changing: %-20s --> %-20s"%(file_name, replace_key, str(val)))
        replace_func( file_name, replace_key, val )
        record[replace_key] = val
        iter_file_parameters( file_name, paras_copy, flag_file, record, rest_files, **kwargs ) # replace the next parameter
    return

def replace_parameter( file_name, para_name, val ):
    with open( file_name, 'r' ) as f:
        content = f.read()

    content_new, Nsub = re.subn( r"^(%s\s+)([^\s]+)"%para_name, r"\g<1>%s"%str(val), content, flags=re.M )
    if Nsub == 0: raise BaseException("ERROR: Cannot find <%s> in <%s>."%(para_name, file_name))

    with open( file_name, 'w' ) as f:
        f.write( content_new )
    return

def replace_parameter_flag( file_name, para_name, val ):
    with open( file_name, 'r' ) as f:
        header = f.readline()

    index = { v:i for i, v in enumerate(header.split()[1:]) } # the first element is assumed to be "#"
    data  = np.loadtxt(file_name)

    if len(index) != data.shape[1]: raise BaseException("ERROR: The number of columns in <%s> does not match the header."%(file_name))

    data[:, index[para_name]] = val

    np.savetxt( file_name, data, fmt="%7g"+"%20.16g"*(len(index)-1), header=header[2:-1] )
    return

def execution( **kwargs ):
    """
    Main execution after iterating all the parameters.
    """
    record_parameters = kwargs["record"]
    print( record_parameters )
    print( "GO GO GAMER GO!" )
    return RETURN_SUCCESS



#====================================================================================================
# Main
#====================================================================================================
if __name__ == "__main__":
    sys.setrecursionlimit( RECURSION_LIMIT ) # reset the recursion depth limit

    # 1. Set up the files and the parameters to be changed
    file_name1   = "Input__Parameter"                                             # file name
    const_paras1 = { "END_T":-1, "END_STEP":50 }                                  # the constant parameters
    iter_paras1  = { "NX0_TOT_X":[16, 32], "NX0_TOT_Y":[16, 32] }                 # the parameters iterated as the given list
    file1        = File( file_name1, const_paras1, iter_paras1, flag_file=False ) # set the `File` class

    file_name2   = "Input__TestProb"
    const_paras2 = { "Const_1":-1, "Const_2":50 }
    iter_paras2  = { "Parameter_1":[0.01, 0.02], "Parameter_2":[16, 32] }
    #file2        = File( file_name2, const_paras2, iter_paras2, flag_file=False )

    # this will change the single column to the same value
    file_name3   = "Input__Flag"
    const_paras3 = { "derefine":0.1 }
    iter_paras3  = { "Soften":[0.01, 0.02], "refine":[16, 32] }
    #file3        = File( file_name3, const_paras3, iter_paras3, flag_file=True )

    # this will change the column to the assigned value
    file_name4   = "Input__Flag_Rho"
    const_paras4 = {}
    iter_paras4  = { "Density":[( 0, 2, 4, 8,10,12,14,16,18,20,22,24 ),
                                ( 1, 3, 5, 7, 9,11,13,15,17,19,21,23 )] }
    file4        = File( file_name4, const_paras4, iter_paras4, flag_file=True )

    files = [file1, file2, file3, file4]        # wrap all the `File` classes.

    # 2. Taking the input arguments
    parser = argparse.ArgumentParser( description = "A script of changing the parameters for Input__* files.",
                                      formatter_class = argparse.RawTextHelpFormatter,
                                      add_help=False)

    parser.add_argument( "-h", "--help",
                         action="help", default=argparse.SUPPRESS,
                         help="Show this help message and exit.\n"
                       )

    parser.add_argument( "-q", "--quite",
                         action="store_true",
                         help="Enable silent mode.\n"
                       )

    parser.add_argument( "-e", "--exe", type=str,
                         default="./gamer",
                         help="The file you want to execute.\n"
                       )

    args = vars( parser.parse_args() )

    # 3. Start changing parameters and running
    iter_files( files, **args )

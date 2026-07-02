#!/bin/python3
"""
A script for changing the parameters of Input__* files.

How to use it:
  1. Set up the files and parameters you want to change in the `Main` section.
     The current script can be used directly for the Plummer test problem.
    a. The file is set to be a `File` class. The `File` takes four inputs:
       file name, constant parameters, variable parameters, and flag file or not.
       Check the `Classes` section for details.
    b. Wrap all the `File` classes to a list called `files`.

    * NOTICE:
      Parameter files:
        The first and second columns specify the name and value of the target parameter, respectively.
      Flag files:
        The first line should always start with `#` and the following columns in the header specify
        the name of each flag threshold.

  2. [Optional] Tailor the `execution()` function for your tests.

  3. [Optional] To pass extra arguments to `execution()` from the command line,
     add your own `parser.add_argument` in the `Main` section.
-----------------------------------------------------------------------------------------------------
For developer:
1. The main concept is to use a recursive function instead of nested for loops so that the code stays clean
   and easy to maintain.
2. We iterate the target files first and then the parameters of each file. At the end of iteration, we call
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
import shutil
import glob



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
        consts    : dict. The constant parameters to be specified.
        paras     : dict with list elements. The parameters to be changed.
        flag_file : bool. Whether or not the target file is a refinement flag file.
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
    Main execution after iterating all parameters.
    """
    cwd               = os.getcwd()
    record_parameters = kwargs["record"]

    # 1. Add parameter change by other parameters
    record_parameters["NX0_TOT_Y"] = record_parameters["NX0_TOT_X"]
    replace_parameter( "Input__Parameter", "NX0_TOT_Y", record_parameters["NX0_TOT_Y"] )
    record_parameters["NX0_TOT_Z"] = record_parameters["NX0_TOT_X"]
    replace_parameter( "Input__Parameter", "NX0_TOT_Z", record_parameters["NX0_TOT_Z"] )

    # 2. Run gamer
    # subprocess.run( ["./gamer > log 2>&1"], shell=True )
    subprocess.run( ["mpirun -map-by ppr:2:socket:pe=8 --report-bindings ./gamer 1>>log 2>&1"], shell=True )
    # subprocess.run( ["mpirun -map-by ppr:4:socket:pe=8 --report-bindings ./gamer 1>>log 2>&1"], shell=True )

    # 3. Analysis: Call python scripts etc.

    # 4. Create a folder named by the parameters to be changed
    par_dir  = "gamer_" + "_".join(record_parameters.keys())

    if not os.path.isdir(par_dir): os.mkdir(par_dir)

    # 5. Create a subfolder named by the current runtime parameters
    sub_dir  = str(record_parameters)
    dest_dir = os.path.join(par_dir, sub_dir)
    if not os.path.isdir(dest_dir): os.mkdir(dest_dir)

    # 6. Move input and output files to the subfolder; delete with os.remove() (files) and shutil.rmtree() (directories) if necessary
    move_files = [ r'*.png', r'Record*', r'Data*', r'Particle_*', r'log']
    copy_files = [ r'Input*', ]

    for f_type in move_files:
        for f in glob.glob(f_type):
            shutil.move(os.path.join(cwd, f), os.path.join(cwd, dest_dir))

    for f_type in copy_files:
        for f in glob.glob(f_type):
            shutil.copy(os.path.join(cwd, f), os.path.join(cwd, dest_dir))
    return RETURN_SUCCESS



#====================================================================================================
# Main
#====================================================================================================
if __name__ == "__main__":
    sys.setrecursionlimit( RECURSION_LIMIT ) # reset the recursion depth limit

    # 1. Set up the files and parameters to be changed
    file_name1   = "Input__Parameter"                                             # file name
    const_paras1 = { "END_T":-1, "END_STEP":5 }                                   # the constant parameters
    iter_paras1  = { "NX0_TOT_X":[64, 128], "MAX_LEVEL":[2, 3] }                  # the parameters iterated as the given list
    file1        = File( file_name1, const_paras1, iter_paras1, flag_file=False ) # set the `File` class

    # this will change the entire column to the same value
    file_name2   = "Input__Flag_NParPatch"
    const_paras2 = {}
    iter_paras2  = { "Number_of_particles_per_patch":[200, 400] }
    file2        = File( file_name2, const_paras2, iter_paras2, flag_file=True )

    # this will change the column to the assigned values
    file_name3   = "Input__Flag_Rho"
    const_paras3 = {}
    iter_paras3  = { "Density":[ tuple( [10**(i-3) for i in range(12)] ),
                                 tuple( [10**(i-4) for i in range(12)] )] }
    file3        = File( file_name3, const_paras3, iter_paras3, flag_file=True )

    files = [file1, file2, file3] # wrap all the `File` classes.

    # 2. Taking the input arguments
    parser = argparse.ArgumentParser( description = "A script for changing the parameters of Input__* files.",
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

    args = vars( parser.parse_args() )

    # 3. Start iterating parameters and running gamer
    iter_files( files, **args )

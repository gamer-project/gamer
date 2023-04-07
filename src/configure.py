"""
A. User Guide:
  This script is for generating the GAMER Makefile. To use this script, you need to know the followings:
  0. Help:
    Use the `-h` or `--help` to get the short help message. Use `-lh` to get detail help message.
  1. Library paths:
    When using library, you need to assign the library path. The path config file can be found under `../configs/`. 
    To setup your own config file, please copy `example.config` and modify it.
  2. Compilation flags: 
    We already have two flag config files under `../configs/`, `intel.make` and `gnu.make`, which are for the 
    intel and gnu compiler respetively. To setup your own config file, please copy `example.make` and modify it.
  3. Run the script: 
    This script can run under Python2 and Python3. To run the script, please run the following command:
      `python configure.py [--your_arguments]`
    After the command finishes, the `Makefile` will be generated.

B. Developer guide:
  This code can be separated to the 5 parts which are `Packages`, `Global variables`, `Classes`, `Functions`, and
  `Main execution`.

  1. Add new simulation option:
    a. Make sure you have added the necessary source files in `Makefile_base`.
    b. Add the python argument reader for the new simulation option in the `Main execution` part of the code.
    c. Add the name convertion dictionary `[python_argument:gamer_argument]` to the `NAME_TABLE`.
       (You can find it in the `Global variable` part of the code.)
    d. Add the new argument to the `sim_opt` parameter of `load_sim()`. (You can find it in the `Functions` 
       part of the code.)
    e. [Optional] Add the error rules in `validation()`. (You can find it in the `Functions` part of the code.)
    f. [Optional] Add the warning rules in `warning()`. (You can find it in the `Functions` part of the code.)
  2. Add a new path:
    a. Add a new line in the Makefile_base in the path section. 
      `NEW_PATH := @@@NEW_PATH@@@`
    b. Add a new line in the path config file. 
      `NEW_PATH    /path/of/new`
  3. Add a new compiler flag type:
    a. Add a new line in the Makefile_base in the flag section. 
      `NEW_FLAG := @@@NEW_FLAG@@@`
    b. Add a new flag key `["NEW_FLAG":""]` in `flags` of `load_compile()`.
    c. Add a new line in the path config file. 
      `NEW_FLAG    -new_flag`
  4. Rules of Makefile_base:
    a. The string will be replaced by this script need to be sandwitched by `@@@`.
  5. Modify the help message:
    We overwrite the original `print_help` since GAMER has too many options. You should follow the followings.
    a. `print_help`: Showing the essential message only. The default option should start with '*' sign.
    b. `print_help_detail`: Should be as clear as possible. 
"""

####################################################################################################
# Packages
####################################################################################################
import argparse
import os
import sys
import re



####################################################################################################
# Global variables
####################################################################################################
GAMER_CONFIG_DIR  = "../configs/"
GAMER_MAKE_BASE   = "Makefile_base"
GAMER_MAKE_OUT    = "Makefile"
GAMER_DESCRIPTION = "Prepare customize Makefile for GAMER.\nThe default value starts with '*' sign.\nTo show the detail help message, please use -lh argument."
GAMER_EPILOG      = "The default complie flags are %sintel.make and %sgnu.make"%(GAMER_CONFIG_DIR, GAMER_CONFIG_DIR)
# The convert name from the python argument to the makefile argument.
NAME_TABLE        = {"model":"MODEL", "passive":"NCOMP_PASSIVE_USER", "flu_scheme":"FLU_SCHEME", 
                     "slope":"LR_SCHEME", "flux":"RSOLVER", "dual":"DUAL_ENERGY", "mhd":"MHD",
                     "cosmic_ray":"COSMIC_RAY", "eos":"EOS", "barotropic":"BAROTROPIC_EOS", 
                     "conserve_mass":"CONSERVE_MASS", "laplacian_four":"LAPLACIAN_4TH", 
                     "self_interaction":"QUARTIC_SELF_INTERACTION", "gravity":"GRAVITY", 
                     "pot_scheme":"POT_SCHEME", "store_pot_ghost":"STORE_POT_GHOST", 
                     "unsplit_gravity":"UNSPLIT_GRAVITY", "comoving":"COMOVING", "particle":"PARTICLE", 
                     "tracer":"TRACER", "store_acc":"STORE_PAR_ACC", "star_formation":"STAR_FORMATION", 
                     "par_attribute":"PAR_NATT_USER", "nlevel":"NLEVEL", "max_patch":"MAX_PATCH", 
                     "patch_size":"PATCH_SIZE", "debug":"GAMER_DEBUG", 
                     "bitwise_reproduce":"BITWISE_REPRODUCIBILITY", "timing":"TIMING", 
                     "timing_solver":"TIMING_SOLVER", "double":"FLOAT8", "laohu":"LAOHU", 
                     "hdf5":"SUPPORT_HDF5", "GSL":"SUPPORT_GSL", "FFTW":"SUPPORT_FFTW", 
                     "LIBYT":"SUPPORT_LIBYT", "LIBYT_patch":"LIBYT_USE_PATCH_GROUP", "RNG":"RANDOM_NUMBER", 
                     "serial_compiler":"SERIAL", "openmp":"OPENMP", "mpi":"LOAD_BALANCE=HILBERT", 
                     "overlap_mpi":"OVERLAP_MPI", "GPU":"GPU", "GPU_arch":"GPU_ARCH"}
class BCOLOR:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'



####################################################################################################
# Classes
####################################################################################################
class ArgumentParser( argparse.ArgumentParser ):
    def __init__(self, *args, **kwargs):
        self.program = { key: kwargs[key] for key in kwargs }
        self.options = []
        super(ArgumentParser, self).__init__(*args, **kwargs)

    def add_argument(self, *args, **kwargs):
        super(ArgumentParser, self).add_argument(*args, **kwargs)
        option = {}
        option["flags"] = [ item for item in args ]
        for key in kwargs:
            option[key] = kwargs[key]
        self.options.append(option)

    def string_align( self, string, indent, width, end_char ):
        """
        end_char : The ending character of a word.
        """
        N          = len(indent)
        sub_indent = N * " "
        if width < N:  raise ValueError("width is smaller than indent length.")
        
        now_n = 0
        new_str = ""
        new_line = False
        for i in range(len(string)):
            if new_line:
                new_str += "\n" + sub_indent
                new_line = False
                now_n = N
                
            if string[i] == "\n": 
                new_line = True
                continue
            
            new_str += string[i]
            now_n += 1
            
            if now_n >= width:
                if string[i] == end_char: new_line = True
        return new_str

    def print_usage(self, *args, **kwargs):
        usage_width  = 100

        if "usage" in self.program:
            print("Usage: %s" % self.program["usage"])
        else:
            usage = []
            for option in self.options:
                for item in option["flags"]:
                    if "choices" in option:
                        if "default" in option:
                            temp = [ "*"+str(opt) if opt == option["default"] else str(opt) for opt in option["choices"] ]
                        else:
                            temp = [ str(opt) for opt in option["choices"] ]
                        usage += [ "[%s {%s}]"%(item, ", ".join(temp)) ]
                        continue
                    
                    if "metavar" in option:
                        if "default" in option:
                            usage += [ "[%s %s *%s]"%(item, option["metavar"], str(option["default"])) ]
                        else:
                            usage += [ "[%s %s]"%(item, option["metavar"]) ]
                        continue
                    
                    if "dest" in option:
                        if "default" in option:
                            usage += [ "[%s %s *%s]"%(item, option["dest"], str(option["default"])) ]
                        else:
                            usage += [ "[%s %s]"%(item, option["dest"].upper()) ]
                        continue
                    
                    if "action" in option:
                        if option["action"] in ["help", "store_const", "store_true", "store_false"]:
                            usage += ["[%s]"%(item)]
                            continue

                    temp = re.sub(r"^(-{1,})", "", item).upper()
                    if "default" in option:
                        usage += [ "[%s %s *%s]"%(item, temp, str(option["default"])) ]
                    else:
                        usage += [ "[%s %s]"%(item, temp) ]
            indent     = "Usage: %s " % os.path.basename(sys.argv[0])
            output = indent + " " + str.join(" ", usage)
            print( self.string_align(output, indent, usage_width, "]") )
        print("")

    def print_help(self, *args, **kwargs):
        # Print usage
        self.print_usage()

        # Print description
        if "description" in self.program:
            print(self.program["description"])
            print("")
        
        # Print epilog
        if "epilog" in self.program:
            print(self.program["epilog"])
            print("")

    def print_help_detail(self):
        # Print usage
        self.print_usage()

        # Print description
        if "description" in self.program:
            print(self.program["description"])
            print("")

        # Print options
        print("Options:")
        option_width = 100
        option_indent = 0
        for option in self.options:
            option["flags2"] = str.join(", ", [ "%s %s" % (item, option["metavar"]) if "metavar" in option else "%s %s" % (item, option["dest"].upper()) if "dest" in option else item for item in option["flags"] ])
            if len(option["flags2"]) > option_indent:
                option_indent = len(option["flags2"])
        
        for option in self.options:
            template = "  %-" + str(option_indent) + "s  "
            indent = template %(option["flags2"])
            output = indent

            if "help" in option: output += option["help"]
            
            if "action" in option:
                if option["action"] == "help":
                    print( self.string_align(output, indent, option_width, " ") )
                    continue
            
            if "choices" in option:
                temp = [ str(opt) for opt in option["choices"] ]
                output += "Choice: [%s] => "%(", ".join(temp))

            if "default" in option:
                output += "Default: %s" % option["default"] if isinstance(option["default"], str) else "Default: %s" % str(option["default"])

            if "action" in option:
                output += "Default: False" if option["action"] == "store_true" else "Default: False"

            print( self.string_align(output, indent, option_width, " ") )
        
        # Print epilog
        if "epilog" in self.program:
            print(self.program["epilog"])
            print("")



####################################################################################################
# Functions 
####################################################################################################
def color_print( string, color ):
    print( color + string + BCOLOR.ENDC )
    return

def load_config( config ):
    print("Using %s as the path config."%(config))
    paths = {}
    with open( config, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        if line[0] == "#": continue
        temp = line.split()
        if len(temp) == 0: continue
        try:
           paths[temp[0]] = temp[1]
        except:
           color_print( "WARNING: %s is not set."%temp[0], BCOLOR.WARNING )
           paths[temp[0]] = ''

    return paths

def load_flag( config ):
    print("Using %s as the flag config."%(config))
    flags = {"CXXFLAG":"", "OPENMPFLAG":"", "LIBFLAG":"", "CUDAFLAG":""}
    with open( config, 'r') as f:
        lines = f.readlines()
   
    for line in lines:
        if line[0] == "#":    continue 
        temp = line.split()
        if len(temp) == 0: continue
        if len(temp) == 1: continue

        for i in range(1, len(temp)):
            if temp[i][0] == "#": break
            flags[temp[0]] += temp[i] + " "

    return flags

def load_sims( **kwargs ):
    sim_opt = {}

    # A. Physics
    # A.1 Module 
    if   kwargs["model"] == "HYDRO":
        sim_opt[NAME_TABLE["model"]]      = kwargs["model"]
        sim_opt[NAME_TABLE["flu_scheme"]] = kwargs["flu_scheme"]
        sim_opt[NAME_TABLE["slope"]]      = kwargs["slope"]
        sim_opt[NAME_TABLE["flux"]]       = kwargs["flux"]
        sim_opt[NAME_TABLE["eos"]]        = "EOS_" + kwargs["eos"]
        sim_opt[NAME_TABLE["passive"]]    = kwargs["passive"]

        if kwargs["mhd"]:         sim_opt[NAME_TABLE["mhd"]]  = kwargs["mhd"]
        if kwargs["dual"] != "":  sim_opt[NAME_TABLE["dual"]] = kwargs["dual"]
        if kwargs["cosmic_ray"] : sim_opt[NAME_TABLE["cosmic_ray"]] = kwargs["cosmic_ray"]
        if kwargs["barotropic"]:  sim_opt[NAME_TABLE["barotropic"]] = kwargs["barotropic"]


    elif kwargs["model"] == "ELBDM":
        sim_opt[NAME_TABLE["model"]]   = kwargs["model"]
        sim_opt[NAME_TABLE["passive"]] = kwargs["passive"]
        if kwargs["conserve_mass"]:    sim_opt[NAME_TABLE["conserve_mass"]]    = kwargs["conserve_mass"]
        if kwargs["laplacian_four"]:   sim_opt[NAME_TABLE["laplacian_four"]]   = kwargs["laplacian_four"]
        if kwargs["self_interaction"]: sim_opt[NAME_TABLE["self_interaction"]] = kwargs["self_interaction"]

    elif kwargs["model"] == "PAR_ONLY":
        sim_opt[NAME_TABLE["model"]] = kwargs["model"]

    # A.2 Gravity
    if kwargs["gravity"]:
        sim_opt[NAME_TABLE["gravity"]]     = kwargs["gravity"]
        sim_opt[NAME_TABLE["pot_scheme"]]  = kwargs["pot_scheme"]
        if kwargs["store_pot_ghost"]: sim_opt[NAME_TABLE["store_pot_ghost"]] = kwargs["store_pot_ghost"]
        if kwargs["unsplit_gravity"]: sim_opt[NAME_TABLE["unsplit_gravity"]] = kwargs["unsplit_gravity"]
        if kwargs["comoving"]       : sim_opt[NAME_TABLE["comoving"]]        = kwargs["comoving"]

    # A.3 Particle
    if kwargs["particle"]:
        sim_opt[NAME_TABLE["particle"]]      = kwargs["particle"]
        sim_opt[NAME_TABLE["par_attribute"]] = kwargs["par_attribute"]
        if kwargs["tracer"]         : sim_opt[NAME_TABLE["tracer"]]         = kwargs["tracer"]
        if kwargs["store_acc"]      : sim_opt[NAME_TABLE["store_acc"]]      = kwargs["store_acc"]
        if kwargs["star_formation"] : sim_opt[NAME_TABLE["star_formation"]] = kwargs["star_formation"]

    # A.4 Grackle
    if kwargs["grackle"] : sim_opt[NAME_TABLE["grackle"]] = kwargs["grackle"]

    # B. miscellaneous options
    sim_opt[NAME_TABLE["nlevel"]]     = kwargs["nlevel"]
    sim_opt[NAME_TABLE["max_patch"]]  = kwargs["max_patch"]
    sim_opt[NAME_TABLE["patch_size"]] = kwargs["patch_size"]
    sim_opt[NAME_TABLE["RNG"]]        = kwargs["RNG"]
    if kwargs["bitwise_reproduce"] : sim_opt[NAME_TABLE["bitwise_reproduce"]] = kwargs["bitwise_reproduce"]
    
    if kwargs["debug"]         : sim_opt[NAME_TABLE["debug"]]         = kwargs["debug"] 
    if kwargs["timing"]        : sim_opt[NAME_TABLE["timing"]]        = kwargs["timing"]
    if kwargs["timing_solver"] : sim_opt[NAME_TABLE["timing_solver"]] = kwargs["timing_solver"]
    if kwargs["double"]        : sim_opt[NAME_TABLE["double"]]        = kwargs["double"] 
    if kwargs["laohu"]         : sim_opt[NAME_TABLE["laohu"]]         = kwargs["laohu"]
    if kwargs["hdf5"]          : sim_opt[NAME_TABLE["hdf5"]]          = kwargs["hdf5"]
    if kwargs["GSL"]           : sim_opt[NAME_TABLE["GSL"]]           = kwargs["GSL"]
    if kwargs["FFTW"]          : sim_opt[NAME_TABLE["FFTW"]]          = kwargs["FFTW"]
    if kwargs["LIBYT"]         : sim_opt[NAME_TABLE["LIBYT"]]         = kwargs["LIBYT"]
    if kwargs["LIBYT_patch"]   : sim_opt[NAME_TABLE["LIBYT_patch"]]   = kwargs["LIBYT_patch"]

    # C. parallel options
    if kwargs["openmp"]: sim_opt[NAME_TABLE["openmp"]] = kwargs["openmp"]
    if kwargs["mpi"]:
        sim_opt[NAME_TABLE["mpi"]] = kwargs["mpi"]
    else:
        sim_opt["SERIAL"] = True  # hard coded the option of serial

    if kwargs["GPU"]:
        sim_opt[NAME_TABLE["GPU"]] = kwargs["GPU"]
        sim_opt[NAME_TABLE["GPU_arch"]] = kwargs["GPU_arch"]

    # D. Setup the simulation option string.
    # NOTE: every -Doption must have trailing space
    opts = ""
    for key, val in sim_opt.items():
        if type(val) == type(True): 
            opts += "-D%s "%(key)
            print("%-20s : %r"%(key, val))
        elif type(val) == type("str"): 
            opts += "-D%s=%s "%(key, val)
            print("%-20s : %s"%(key, val))
        elif type(val) == type(0):
            opts += "-D%s=%d "%(key, val)
            print("%-20s : %d"%(key, val))
        else:
            raise TypeError("Unknown type to add the simulation options.")

    return {"SIMU_OPTION":opts}

def load_compile( paths, flags, kwargs ):
    com_opt = {}

    # 1. complier. If mpi, overwrite the compiler to mpi.
    com_opt["CXX"] = paths["MPI_PATH"] + "/bin/mpicxx" if kwargs["mpi"] else kwargs["serial_compiler"]

    # 2. compiler flags
    if kwargs["serial_compiler"] == "icpc" and not kwargs["openmp"]:
        flags["CXXFLAG"] += "-Wno-unknown-pragmas -diag-disable 3180 "

    elif kwargs["serial_compiler"] == "g++" and not kwargs["openmp"]:
        flags["CXXFLAG"] += "-Wno-unknown-pragmas "
    elif not kwargs["openmp"]:
        flags["CXXFLAG"] += "-Wno-unknown-pragmas "

    # 3. openmp flags
    if not kwargs["openmp"]: flags["OPENMPFLAG"] = ""

    # 4. write to complie option
    for key, val in flags.items():
        com_opt[key] = val

    return com_opt

def validation( paths, **kwargs ):
    """
    validation then store the parameter
    """
    success = True
    
    # 0. Makefile 
    if not os.path.isfile( GAMER_MAKE_BASE ):
        color_print("ERROR: %s does not exist."%(GAMER_MAKE_BASE), BCOLOR.FAIL)
        success = False
    
    sim_opt = {}
    
    # A. Physics
    # A.1 Module 
    if kwargs["model"] == "HYDRO":
        if kwargs["mhd"]: 
            if kwargs["flu_scheme"] not in ["MHM_RP", "CTU"]:
                color_print("MHD only supports MHM_RP and CTU.", BCOLOR.FAIL)
                success = False
            if kwargs["flux"] not in ["EXACT", "ROE", "HLLE", "HLLD"]:
                color_print("MHD only supports EXACT, ROE, HLLE, and HLLD Riemann solver.", BCOLOR.FAIL)
                success = False
        else:
            if kwargs["flux"] not in ["EXACT", "ROE", "HLLE", "HLLC"]:
                color_print("Pure hydro only supports EXACT, ROE, HLLE, and HLLC Riemann solver.", BCOLOR.FAIL)
                success = False

        if kwargs["flu_scheme"] == "RTVD":
            if kwargs["passive"] != 0: 
                color_print("Passive scalar is not supported for RTVD.", BCOLOR.FAIL)
                success = False
            if kwargs["unsplit_gravity"]: 
                color_print("Unsplit gravity is not supported for RTVD.", BCOLOR.FAIL)
                success = False
            if kwargs["dual"] != "": 
                color_print("Dual energy is not supported for RTVD.", BCOLOR.FAIL)
                success = False

        if kwargs["dual"] not in ["", "ENPY"]: 
            color_print("This dual energy form is not supported yet.", BCOLOR.FAIL)
            success = False
        
        if kwargs["eos"] != "GAMMA":
            if kwargs["dual"] == "ENPY":
                color_print("ENPY dual energy is only supported for --eos=GAMMA.", BCOLOR.FAIL)
                success = False
            if kwargs["flux"] != "ROE" or kwargs["flux"] != "EXACT":
                color_print("Only ROE and EXACT Riemann solver are supported for --eos!=GAMMA.", BCOLOR.FAIL)
                success = False
            if kwargs["flu_scheme"] == "RTVD" or kwargs["flu_scheme"] == "CTU":
                color_print("RTVD and CTU are only supported for --eos=GAMMA.", BCOLOR.FAIL)
                success = False
            if kwargs["comoving"]:
                color_print("--comoving is only supported for --eos=GAMMA.", BCOLOR.FAIL)
                success = False

        if kwargs["eos"] == "ISOTHERMAL" and not kwargs["barotropic"]:
            color_print("--barotropic must be enabled for --eos=ISOTHERMAL.", BCOLOR.FAIL)
            success = False
            
        if kwargs["cosmic_ray"] : 
            color_print("--cosmic_ray is not supported yet.", BCOLOR.FAIL)
            success = False
            if kwargs["dual"] == "ENPY":
                color_print("ENPY dual energy is not supported for --cosmic_ray.", BCOLOR.FAIL)
                success = False
        
        if kwargs["barotropic"]:
            if kwargs["eos"] not in ["ISOTHERMAL", "TABULAR", "USER"]:
                color_print("--barotropic is only supported for ISOTHERMAL, TABULAR, and USER.", BCOLOR.FAIL)
                success = False

    elif kwargs["model"] == "ELBDM":
        if kwargs["self_interaction"]: 
            if not kwargs["gravity"]:
                color_print("--gravity must be enabled when --self_interaction.", BCOLOR.FAIL)
                success = False
            if kwargs["comoving"]:
                color_print("--comoving is not supported with --self_interaction.", BCOLOR.FAIL)
                success = False
        if kwargs["passive"] != 0:
            color_print("Not supported yet and can only be used as auxiliary fields.", BCOLOR.FAIL)
            success = False

    elif kwargs["model"] == "PAR_ONLY":
        color_print("--model PAR_ONLY is not supported yet.", BCOLOR.FAIL)
        success = False
    else:
        color_print("Unrecognized model: %s. Please add to the model choices."%kwargs["model"], BCOLOR.FAIL)
        success = False

    # A.2 Gravity
    if kwargs["gravity"]:
        if not kwargs["FFTW"]:
            color_print("--FFTW must be enable with --gravity.", BCOLOR.FAIL)
            success = False
        if kwargs["unsplit_gravity"] and kwargs["model"] != "HYDRO":
            color_print("--unsplit_gravity is only supported for --model=HYDRO.", BCOLOR.FAIL)
            success = False
     
    # A.3 Particle
    if kwargs["particle"]:
        if kwargs["star_formation"] and kwargs["store_par_acc"]:
            if not kwargs["store_pot_ghost"]:
                color_print("--store_pot_ghost must be enabled when --star_formation and --store_par_acc are enabled.", BCOLOR.FAIL)
                success = False
        if not kwargs["gravity"] and not kwargs["tracer"]:
            color_print("One of --gravity or --tracer must be enabled for --particle.", BCOLOR.FAIL)
            success = False

    # A.4 Grackle
    if kwargs["grackle"] and kwargs["eos"] != "GAMMA":
        color_print("--grackle must work with --eos=GAMMA.", BCOLOR.FAIL)
        success = False

    # B. miscellaneous options
    if kwargs["nlevel"] < 1:
        color_print("--nlevel should be greater than zero.", BCOLOR.FAIL)
        success = False

    if kwargs["max_patch"] < 1:
        color_print("--max_patch should be greater than zero.", BCOLOR.FAIL)
        success = False

    if kwargs["patch_size"]%2 != 0 or kwargs["patch_size"] < 8:
        color_print("--patch_size should be even number and greater than 6.", BCOLOR.FAIL)
        success = False

    if not kwargs["timing"] and kwargs["timing_solver"]:
        color_print("--timing_solver must enable --timing.", BCOLOR.FAIL)
        success = False

    if kwargs["LIBYT_patch"] and not kwargs["LIBYT"]:
        color_print("--LIBYT_patch must enable --LIBYT.", BCOLOR.FAIL)
        success = False

    if kwargs["overlap_mpi"]:
        color_print("--overlap_mpi is not supported yet.", BCOLOR.FAIL)
        success = False
        if not kwargs["mpi"]:
            color_print("--overlap_mpi must enable --mpi.", BCOLOR.FAIL)
            success = False
    
    if not success: raise BaseException(BCOLOR.FAIL+"The above validations are fail."+BCOLOR.ENDC)
    return


def warning( paths, **kwargs ):
    # 0. Makefile 
    if os.path.isfile( GAMER_MAKE_OUT ):
        color_print("Warning: %s already exist and will be overwriten."%(GAMER_MAKE_OUT), BCOLOR.WARNING)

    # 1. serial config not match
    if kwargs["serial_compiler"] == "icpc" and kwargs["flags"] == "gnu":
        color_print("Warning: The compiler does not match to the default flag config.", BCOLOR.WARNING)
    
    if kwargs["serial_compiler"] == "g++" and kwargs["flags"] == "intel":
        color_print("Warning: The compiler does not match to the default flag config.", BCOLOR.WARNING)


    # 2. mpi
    if kwargs["mpi"]:
        color_print("Warning: MPI is set, Serial mode force to be disabled.", BCOLOR.WARNING)
    
    
    # 3. RTVD
    if kwargs["flu_scheme"] == "RTVD":
        color_print("Warning: Data reconstruction is useless for RTVD.", BCOLOR.WARNING)
        color_print("Warning: Riemann solver is useless for RTVD.", BCOLOR.WARNING)
    
    # 4. RNG sys
    if kwargs["RNG"] == "RNG_GNU_EXT":
        color_print("Warning: %s is not supported on some macOS."%kwargs["RNG"], BCOLOR.WARNING)
    
    if kwargs["RNG"] == "RNG_CPP":
        color_print("Warning: You may need to add -std=c++11 to CXXFLAG in flag config.", BCOLOR.WARNING)
    
    # 5. Path
    if kwargs["GPU"]:
        if "CUDA_PATH" not in paths:
            color_print("CUDA_PATH is not given with --GPU.", BCOLOR.WARNING)
            paths["CUDA_PATH"] = ""
        elif paths["CUDA_PATH"] == "":
            color_print("CUDA_PATH is not given with --GPU.", BCOLOR.WARNING)
    
    if kwargs["FFTW"]:
        if "FFTW_PATH" not in paths:
            color_print("FFTW_PATH is not given with --FFTW.", BCOLOR.WARNING)
            paths["FFTW_PATH"] = ""
        elif paths["FFTW_PATH"] == "":
            color_print("FFTW_PATH is not given with --FFTW.", BCOLOR.WARNING)
    
    if kwargs["mpi"]:
        if "MPI_PATH" not in paths:
            color_print("MPI_PATH is not given with --mpi.", BCOLOR.WARNING)
            paths["MPI_PATH"] = ""
        elif paths["MPI_PATH"] == "":
            color_print("MPI_PATH is not given with --mpi.", BCOLOR.WARNING)
    
    if kwargs["hdf5"]:
        if "HDF5_PATH" not in paths:
            color_print("HDF5_PATH is not given with --hdf5.", BCOLOR.WARNING)
            paths["HDF5_PATH"] = ""
        elif paths["HDF5_PATH"] == "":
            color_print("HDF5_PATH is not given with --hdf5.", BCOLOR.WARNING)
    
    if kwargs["grackle"]:
        if "GRACKLE_PATH" not in paths:
            color_print("GRACKLE_PATH is not given with --grackle.", BCOLOR.WARNING)
            paths["GRACKLE_PATH"] = ""
        elif paths["GRACKLE_PATH"] == "":
            color_print("GRACKLE_PATH is not given with --grackle.", BCOLOR.WARNING)
    
    if kwargs["GSL"]:
        if "GSL_PATH" not in paths:
            color_print("GSL_PATH is not given with --GSL.", BCOLOR.WARNING)
            paths["GSL_PATH"] = ""
        elif paths["GSL_PATH"] == "":
            color_print("GSL_PATH is not given with --GSL.", BCOLOR.WARNING)
    
    if kwargs["LIBYT"]:
        if "LIBYT_PATH" not in paths:
            color_print("LIBYT_PATH is not given with --LIBYT.", BCOLOR.WARNING)
            paths["LIBYT_PATH"] = ""
        elif paths["LIBYT_PATH"] == "":
            color_print("LIBYT_PATH is not given with --LIBYT.", BCOLOR.WARNING)
    return



####################################################################################################
# Main execution
####################################################################################################
# 1. Load the input arguments
parser = ArgumentParser( description = GAMER_DESCRIPTION, 
                         formatter_class = argparse.RawTextHelpFormatter,
                         epilog = GAMER_EPILOG,
                         add_help=False)

parser.add_argument( "-h", "--help",
                     action="help", default=argparse.SUPPRESS,
                     help="Show this help message and exit.\n"
                   )

# detail help message
parser.add_argument( "-lh",
                     action="store_true",
                     help="Show this help message in detail and exit.\n"
                   )

# cluster and flags setup
parser.add_argument( "--cluster", type=str, metavar="NAME",
                     default="eureka",
                     help="Select the cluster. \nChoice: [eureka, YOUR_CLUSTER_NAME] => "
                   )

parser.add_argument( "--flags", type=str, metavar="NAME",
                     default="intel",
                     help="Compiler flags. \nChoice: [intel, gnu, YOUR_FLAG_NAME] => "
                   )

# A. physical models and options of diffierent physical models 
parser.add_argument( "--model", type=str, metavar="MODEL",
                     default="HYDRO", choices=["HYDRO", "ELBDM", "PAR_ONLY"],
                     help="Select the physical model.\n"
                   )

parser.add_argument( "--passive", type=int, metavar="NUMBER",
                     default=0,
                     help="Set the number of passive scalars.\n"
                   )

# A.1 Hydro options
parser.add_argument( "--flu_scheme", type=str, metavar="SCHEME",
                     default="CTU", choices=["RTVD", "MHM", "MHM_RP", "CTU"],
                     help="Select the fluid solver for HYDRO model.\n"
                   )

parser.add_argument( "--slope", type=str, metavar="TYPE",
                     default="PPM", choices=["PLM", "PPM"],
                     help="Select the spatial data reconstruction method.\n"
                   )

parser.add_argument( "--flux", type=str, metavar="TYPE",
                     default="ROE", choices=["EXACT", "ROE", "HLLE", "HLLC", "HLLD"],
                     help="Select the Riemann solver.\n"
                   )

parser.add_argument( "--dual", type=str, metavar="TYPE",
                     default="", choices=["", "ENPY", "EINT"],
                     help="Select the dual-energy formalism.\n"
                   )

parser.add_argument( "--mhd",
                     action="store_true",
                     help="Enable magnetohydrodynamic.\n"
                   )

parser.add_argument( "--cosmic_ray",
                     action="store_true",
                     help="Enable cosmic rays.\n"
                   )

parser.add_argument( "--eos", type=str, metavar="TYPE",
                     default="GAMMA", choices=["GAMMA", "ISOTHERMAL", "NUCLEAR", "TABULAR", "USER"],
                     help="Select the equation of state.\n"
                   )

parser.add_argument( "--barotropic",
                     action="store_true",
                     help="Whether or not the --eos set is barotropic.\n"
                   )

# A.2 ELBDM scheme
parser.add_argument( "--conserve_mass",
                     action="store_true",
                     help="Enforce the mass conservation.\n"
                   )

parser.add_argument( "--laplacian_four",
                     action="store_true",
                     help="Enable the fourth order of Laplacian.\n"
                   )

parser.add_argument( "--self_interaction",
                     action="store_true",
                     help="Including the quartic self-interaction potential.\n"
                   )

# A.3 gravity
parser.add_argument( "--gravity",
                     action="store_true",
                     help="Enable gravity.\n"
                   )

parser.add_argument( "--pot_scheme", type=str, metavar="SCHEME",
                     default="SOR", choices=["SOR", "MG"],
                     help="Select the Poisson solver.\n"
                   )

parser.add_argument( "--store_pot_ghost",
                     action="store_true",
                     help="Store the potential ghost-zone for each patch on each side.\n"
                   )

parser.add_argument( "--unsplit_gravity",
                     action="store_true",
                     help="Use unsplitting method to couple gravity to the target model.\n"
                   )

parser.add_argument( "--comoving",
                     action="store_true",
                     help="Comoving frame for cosmological simulation.\n"
                   )

# A.4 particle
parser.add_argument( "--particle",
                     action="store_true",
                     help="Enable particle.\n"
                   )
parser.add_argument( "--tracer",
                     action="store_true",
                     help="Enable tracer particles.\n"
                   )

parser.add_argument( "--store_acc",
                     action="store_true",
                     help="Store particle acceleration.\n"
                   )

parser.add_argument( "--star_formation",
                     action="store_true",
                     help="Allow creating new particles after initialization.\n"
                   )

parser.add_argument( "--par_attribute", type=int, metavar="NUMBER",
                     default=0,
                     help="Set the number of user defined particle attributes.\n"
                   )

# A.5 grackle
parser.add_argument( "--grackle",
                     action="store_true",
                     help="Enable Grackle, a chemistry and radiative cooling library.\n"
                   )

# B. miscellaneous options
parser.add_argument( "--nlevel", type=int,
                     default=10,
                     help="Set the maximum level of AMR.\n"
                   )

parser.add_argument( "--max_patch", type=int,
                     default=100000,
                     help="Set the maximum patchs on each level of AMR.\n"
                   )

parser.add_argument( "--patch_size", type=int,
                     default=8,
                     help="Set size of each direction of a single patch.\n"
                   )

parser.add_argument( "--debug",
                     action="store_true",
                     help="Enable debug mode.\n"
                   )

parser.add_argument( "--bitwise_reproduce",
                     action="store_true",
                     help="Enable bitwise reproducibility.\n"
                   )

parser.add_argument( "--timing",
                     action="store_true",
                     help="Enable to measure timing.\n"
                   )

parser.add_argument( "--timing_solver",
                     action="store_true",
                     help="Enable to measure GPU time.\n"
                   )

parser.add_argument( "--double",
                     action="store_true",
                     help="Enable double precision.\n"
                   )

parser.add_argument( "--laohu",
                     action="store_true",
                     help="Work on the NAOC Laohu GPU cluster.\n"
                   )

parser.add_argument( "--hdf5",
                     action="store_true",
                     help="Support HDF5 format.\n"
                   )

parser.add_argument( "--GSL",
                     action="store_true",
                     help="Support GNU scientific library.\n"
                   )

parser.add_argument( "--FFTW",
                     action="store_true",
                     help="Support FFTW library.\n"
                   )

parser.add_argument( "--LIBYT",
                     action="store_true",
                     help="Support yt inline analysis.\n"
                   )

parser.add_argument( "--LIBYT_patch",
                     action="store_true",
                     help="Use patch group as the unit in libyt. Note that this will speed up inline-analysis but increase memory consumption.\n"
                   )

parser.add_argument( "--RNG", type=str, metavar="TYPE",
                     default="RNG_GNU_EXT",
                     choices=["RNG_GNU_EXT", "RNG_CPP11"],
                     help="Select the random number generator.\n"
                   )

# C. parallelization and flags
parser.add_argument( "--serial_compiler", type=str, metavar="COMPILER",
                     default="icpc",
                     help="Serial compiler type.\n"
                   )
parser.add_argument( "--openmp",
                     action="store_true",
                     help="Enable OpenMP parallization.\n"
                   )
 
parser.add_argument( "--mpi",
                     action="store_true",
                     help="Enable MPI parallization.\n"
                   )

parser.add_argument( "--overlap_mpi",
                     action="store_true",
                     help="Overlap MPI communication with computation.\n"
                   )

parser.add_argument( "--GPU",
                     action="store_true",
                     help="Enable GPU.\n"
                   )

parser.add_argument( "--GPU_arch", type=str, metavar="TYPE",
                     default="TURING", choices=["FERMI", "KEPLER", "MAXWELL", "PASCAL", "VOLTA", "TURING", "AMPERE"],
                     help="Select the archtecture of GPU.\n"
                   )


args = vars( parser.parse_args() )

# 1.1 print out the detail help message
if args["lh"]:
    parser.print_help_detail()
    exit()


#------------------------------------------------------------
# 2. Prepare the makiefile args
# 2.1 load the cluster steup
paths = load_config( "%s%s.config"%(GAMER_CONFIG_DIR, args["cluster"]) )
flags = load_flag("%s%s.make"%(GAMER_CONFIG_DIR, args["flags"]))

# 2.2 Check if the argument are valid
validation( paths, **args )

warning( paths, **args )

# 2.3 add the SIMU_OPTION
print("")
print("========================================")
print("GAMER has the following setting.")
print("----------------------------------------")
sims = load_sims( **args )

# 2.4 setup compiler
compiles = load_compile( paths, flags, args )


#------------------------------------------------------------
# 3. Create Makefile
# 3.1 Read
with open( GAMER_MAKE_BASE, "r" ) as make_base:
    makefile = make_base.read()

# 3.2 Replace
print("----------------------------------------")
for key, val in paths.items():
    print("%-15s : %s"%(key, val))
    makefile = re.sub(r"@@@%s@@@"%(key), val, makefile)

for key, val in sims.items():
    makefile = re.sub(r"@@@%s@@@"%(key), val, makefile)

print("----------------------------------------")
for key, val in compiles.items():
    print("%-10s : %s"%(key, val))
    makefile = re.sub(r"@@@%s@@@"%(key), val, makefile)

# 3.3 Write
with open( GAMER_MAKE_OUT, "w") as make_out:
    make_out.write( makefile )


print("========================================")
print("%s is created."%GAMER_MAKE_OUT)
print("========================================")

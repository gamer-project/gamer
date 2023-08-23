"""
User and developer guides of this script are provided in the following link.

   https://github.com/gamer-project/gamer/wiki/Installation%3A-Configure.py

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
NONE_STR   = "OFF"
PYTHON_VER = [sys.version_info.major, sys.version_info.minor]

GAMER_CONFIG_DIR  = "../configs"
GAMER_MAKE_BASE   = "Makefile_base"
GAMER_MAKE_OUT    = "Makefile"
GAMER_DESCRIPTION = "Prepare a customized Makefile for GAMER.\nDefault values are marked by '*'.\nUse -lh to show a detailed help message.\n"
GAMER_EPILOG      = "2023 Computational Astrophysics Lab, NTU. All rights reserved.\n"



####################################################################################################
# Classes
####################################################################################################
class BCOLOR:
    HEADER    = '\033[95m'
    OKBLUE    = '\033[94m'
    OKCYAN    = '\033[96m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'

class ArgumentParser( argparse.ArgumentParser ):
    def __init__(self, *args, **kwargs):
        self.program     = { key: kwargs[key] for key in kwargs }
        self.options     = []
        self.depends     = {}
        self.constraints = {}
        self.gamer_names = {}
        # This feature is only supported for version >= 3.5.
        if PYTHON_VER[0] == 2 or PYTHON_VER[1] < 5: kwargs.pop("allow_abbrev")
        super(ArgumentParser, self).__init__(*args, **kwargs)

    def add_argument(self, *args, **kwargs):
        if "depend" in kwargs:
            key = args[0].replace("-", "")
            self.depends[key] = kwargs.pop("depend")
        if "constraint" in kwargs:
            key = args[0].replace("-", "")
            self.constraints[key] = kwargs.pop("constraint")
        if "gamer_name" in kwargs:
            key = args[0].replace("-", "")
            self.gamer_names[key] = kwargs.pop("gamer_name")

        super(ArgumentParser, self).add_argument(*args, **kwargs)
        option = {}
        option["flags"] = [ item for item in args ]
        for key in kwargs:
            option[key] = kwargs[key]
        self.options.append(option)

    def parse_args(self, args=None, namespace=None):
        args, argv = self.parse_known_args(args, namespace)
        msg = "\n"
        close_dist = 2

        for arg in argv:
            if arg[0] != "-":
                msg += 'Unrecognized positional argument: %s\n'%(arg)
                continue
            arg = arg.split("=")[0]     # separate the assigned value.
            min_dist = 100000
            pos_key = ""
            for key in self.gamer_names:
                dist = distance( arg, "--"+key )
                if dist >= min_dist: continue
                min_dist = dist
                pos_key = "--"+key
            msg += 'Unrecognized argument: %s'%(arg)
            msg += ', do you mean: %s ?\n'%(pos_key) if min_dist <= close_dist else "\n"

        if len(argv) != 0: self.error( msg )
        return args, self.gamer_names, self.depends, self.constraints

    def string_align( self, string, indent, width, end_char ):
        """
        end_char : The ending character of a word.
        """
        N          = len(indent)
        sub_indent = N * " "
        if width < N:  raise ValueError("Width is smaller than indent length.")

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
                            usage += [ "[%s %s *%s]"%(item, option["metavar"], "Depend" if option["default"] == None else str(option["default"])) ]
                        else:
                            usage += [ "[%s %s]"%(item, option["metavar"]) ]
                        continue

                    if "dest" in option:
                        if "default" in option:
                            usage += [ "[%s %s *%s]"%(item, option["dest"], "Depend" if option["default"] == None else str(option["default"])) ]
                        else:
                            usage += [ "[%s %s]"%(item, option["dest"].upper()) ]
                        continue

                    if "action" in option:
                        if option["action"] in ["help", "store_const", "store_true", "store_false"]:
                            usage += [ "[%s]"%(item) ]
                            continue

                    temp = re.sub(r"^(-{1,})", "", item).upper()
                    if "default" in option:
                        usage += [ "[%s %s *%s]"%(item, temp, "Depend" if option["default"] == None else str(option["default"])) ]
                    else:
                        usage += [ "[%s %s]"%(item, temp) ]
            indent = "Usage: %s " % os.path.basename(sys.argv[0])
            output = indent + " " + str.join(" ", usage)
            print( self.string_align(output, indent, usage_width, "]") )
        print("")

    def print_help(self, *args, **kwargs):
        # Print usage
        self.print_usage()

        # Print description
        if "description" in self.program: print(self.program["description"])

        # Print epilog
        if "epilog" in self.program: print(self.program["epilog"])

    def print_help_detail(self):
        # Print usage
        self.print_usage()

        # Print description
        if "description" in self.program: print(self.program["description"])

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
                output += "Default: %s" %("Depend" if option["default"] == None else str(option["default"]))

            if "action" in option:
                output += "Default: False" if option["action"] == "store_true" else "Default: False"

            print( self.string_align(output, indent, option_width, " ") )

        # Print epilog
        if "epilog" in self.program: print(self.program["epilog"])



####################################################################################################
# Functions
####################################################################################################
def str2bool(v):
    if isinstance(v, bool): return v
    if v.lower() == "true":
        return True
    elif v.lower() == "false":
        return False
    else:
        raise TypeError("Can not convert <%s> to boolean."%(v))
    return

def color_print( string, color ):
    print( color + string + BCOLOR.ENDC )
    return

def add_option( opt_str, name, val ):
    # NOTE: Every -Doption must have a trailing space.
    if type(val) == type(True):
        if val: opt_str += "-D%s "%(name)
        print("%-25s : %r"%(name, val))
    elif type(val) == type("str"):
        if val != NONE_STR:
            opt_str += "-D%s=%s "%(name, val)
            print("%-25s : %s"%(name, val))
    elif type(val) == type(0):
        opt_str += "-D%s=%d "%(name, val)
        print("%-25s : %d"%(name, val))
    elif type(val) == type(0.):
        opt_str += "-D%s=%f "%(name, val)
        print("%-25s : %f"%(name, val))
    else:
        raise TypeError("The simulation option <%s> has an unknown type <%s>."%(name, str(type(val))))

    return opt_str

def distance( s1, s2 ):
    """
    Calculate the distance between two strings.
    See: https://en.wikipedia.org/wiki/Damerau-Levenshtein_distance
         https://www.geeksforgeeks.org/damerau-levenshtein-distance/
    """
    matrix = [ [ 0 for i in range(len(s2)+1) ] for j in range(len(s1)+1) ]

    for i in range(len(s1)+1):
        matrix[i][0] = i
    for j in range(len(s2)+1):
        matrix[0][j] = j

    for i in range(len(s1)+1):
        for j in range(1, len(s2)+1):
            if s1[i-1] == s2[j-1]:
                matrix[i][j] = matrix[i-1][j-1]
            else:
                matrix[i][j] = 1 + min(matrix[i-1][j], matrix[i][j-1], matrix[i-1][j-1])

    return matrix[len(s1)][len(s2)]

def load_arguments():
    parser = ArgumentParser( description = GAMER_DESCRIPTION,
                             formatter_class = argparse.RawTextHelpFormatter,
                             epilog = GAMER_EPILOG,
                             allow_abbrev=False,
                             add_help=False)

    parser.add_argument( "-h", "--help",
                         action="help", default=argparse.SUPPRESS,
                         help="Show this help message and exit.\n"
                       )

    # detailed help message
    parser.add_argument( "-lh",
                         action="store_true",
                         help="Show this help message in detail and exit.\n"
                       )

    # machine config setup
    parser.add_argument( "--machine", type=str, metavar="MACHINE",
                         default="eureka_intel",
                         help="Select the MACHINE.config file under ../configs directory. \nChoice: [eureka_intel, YOUR_MACHINE_NAME] => "
                       )

    # A. physical models and options of diffierent physical models
    parser.add_argument( "--model", type=str, metavar="TYPE", gamer_name="MODEL",
                         default="HYDRO", choices=["HYDRO", "ELBDM", "PAR_ONLY"],
                         help="The physical model (HYDRO: hydrodynamics/magnetohydrodynamics, ELBDM: wave dark matter, PAR_ONLY: partivle-only). Must be set in any cases. PAR_ONLY is not supported yet.\n"
                       )

    parser.add_argument( "--passive", type=int, metavar="INTEGER", gamer_name="NCOMP_PASSIVE_USER",
                         default=0,
                         depend={"model":["HYDRO", "ELBDM"]},
                         help="Set the number of user-defined passively advected scalars. Useless for RTVD. <--model=ELBDM> doesn't support passive scalars and only regards them as auxiliary fields.\n"
                       )

    # A.1 Hydro options
    parser.add_argument( "--flu_scheme", type=str, metavar="TYPE", gamer_name="FLU_SCHEME",
                         default="CTU", choices=["RTVD", "MHM", "MHM_RP", "CTU"],
                         depend={"model":"HYDRO"},
                         constraint={ "RTVD":{"unsplit_gravity":False, "passive":0, "dual":NONE_STR, "eos":"GAMMA"},
                                      "CTU":{"eos":"GAMMA"} },
                         help="The hydrodynamic/MHD integrator. MHD only supports MHM_RP and CTU.\n"
                       )

    parser.add_argument( "--slope", type=str, metavar="TYPE", gamer_name="LR_SCHEME",
                         default="PPM", choices=["PLM", "PPM"],
                         depend={"model":"HYDRO"},
                         help="The spatial data reconstruction method (PLM: piecewise-linear, PPM: piecewise-parabolic). Useless for <--flu_scheme=RTVD>.\n"
                       )

    parser.add_argument( "--flux", type=str, metavar="TYPE", gamer_name="RSOLVER",
                         default="ROE", choices=["EXACT", "ROE", "HLLE", "HLLC", "HLLD"],
                         depend={"model":"HYDRO"},
                         constraint={ "ROE":{"eos":"GAMMA"},
                                      "EXACT":{"eos":"GAMMA"} },
                         help="The Riemann solver. Pure hydro: EXACT/ROE/HLLE/HLLC^, MHD: ROE/HLLE/HLLD^ (^ indicates the recommended solvers). Useless for RTVD.\n"
                       )

    parser.add_argument( "--dual", type=str, metavar="TYPE", gamer_name="DUAL_ENERGY",
                         default=NONE_STR, choices=[NONE_STR, "DE_ENPY", "DE_EINT"],
                         depend={"model":"HYDRO"},
                         constraint={ "DE_ENPY":{"eos":"GAMMA"} },
                         help="The dual-energy formalism (DE_ENPY: entropy, DE_EINT: internal energy). DE_EINT is not supported yet. Useless for RTVD.\n"
                       )

    parser.add_argument( "--mhd", type=str2bool, metavar="BOOLEAN", gamer_name="MHD",
                         default=False,
                         depend={"model":"HYDRO"},
                         constraint={ True:{"flu_scheme":["MHM_RP", "CTU"], "flux":["ROE", "HLLE", "HLLD"]},
                                     False:{"flux":["EXACT", "ROE", "HLLE", "HLLC"]} },
                         help="Magnetohydrodynamics.\n"
                       )

    parser.add_argument( "--cosmic_ray", type=str2bool, metavar="BOOLEAN", gamer_name="COSMIC_RAY",
                         default=False,
                         depend={"model":"HYDRO"},
                         constraint={ True:{"dual":[NONE_STR]} },
                         help="Cosmic rays (not supported yet).\n"
                       )

    parser.add_argument( "--eos", type=str, metavar="TYPE", gamer_name="EOS",
                         default="GAMMA", choices=["GAMMA", "ISOTHERMAL", "NUCLEAR", "TABULAR", "USER"],
                         depend={"model":"HYDRO"},
                         constraint={ "ISOTHERMAL":{"barotropic":True} },
                         help="Equation of state. Must be set when <--model=HYDRO>. Must enable <--barotropic> for ISOTHERMAL.\n"
                       )

    parser.add_argument( "--barotropic", type=str2bool, metavar="BOOLEAN", gamer_name="BAROTROPIC_EOS",
                         default=False,
                         depend={"model":"HYDRO"},
                         constraint={ True:{"eos":["ISOTHERMAL", "TABULAR", "USER"]} },
                         help="Whether or not the equation of state set by <--eos> is barotropic. Mandatory for <--eos=ISOTHEMAL>. Optional for <--eos=TABULAR> and <--eos=USER>.\n"
                       )

    # A.2 ELBDM scheme
    parser.add_argument( "--conserve_mass", type=str2bool, metavar="BOOLEAN", gamer_name="CONSERVE_MASS",
                         default=True,
                         depend={"model":"ELBDM"},
                         help="Enforce the mass conservation for <--model=ELBDM>.\n"
                       )

    parser.add_argument( "--laplacian_four", type=str2bool, metavar="BOOLEAN", gamer_name="LAPLACIAN_4TH",
                         default=True,
                         depend={"model":"ELBDM"},
                         help="Enable the fourth-order Laplacian for <--model=ELBDM>.\n"
                       )

    parser.add_argument( "--self_interaction", type=str2bool, metavar="BOOLEAN", gamer_name="QUARTIC_SELF_INTERACTION",
                         default=False,
                         depend={"model":"ELBDM"},
                         constraint={ True:{"gravity":True, "comoving":False} },
                         help="Include the quartic self-interaction potential for <--model=ELBDM>. Must enable <--gravity>. Does not support <--comoving>.\n"
                       )

    # A.3 gravity
    parser.add_argument( "--gravity", type=str2bool, metavar="BOOLEAN", gamer_name="GRAVITY",
                         default=False,
                         constraint={ True:{"fftw":["FFTW2", "FFTW3"]} },
                         help="Enable gravity. Must enable <--fftw>.\n"
                       )

    parser.add_argument( "--pot_scheme", type=str, metavar="TYPE", gamer_name="POT_SCHEME",
                         default="SOR", choices=["SOR", "MG"],
                         depend={"gravity":True},
                         help="Select the Poisson solver. SOR: successive-overrelaxation (recommended), MG: multigrid. Must be set when <--gravity> is enabled.\n"
                       )

    parser.add_argument( "--store_pot_ghost", type=str2bool, metavar="BOOLEAN", gamer_name="STORE_POT_GHOST",
                         default=True,
                         depend={"gravity":True},
                         help="Store GRA_GHOST_SIZE ghost-zone potential for each patch on each side. Recommended when PARTICLE is enabled for improving accuaracy for particles around the patch boundaries. Must be enabled for <--star_formation> + <--store_par_acc>.\n"
                       )

    parser.add_argument( "--unsplit_gravity", type=str2bool, metavar="BOOLEAN", gamer_name="UNSPLIT_GRAVITY",
                         default=None,
                         depend={"gravity":True},
                         constraint={ True:{"model":"HYDRO"} },
                         help="Use unsplitting method to couple gravity to the target model (recommended). Supported only for <--model=HYDRO>.\n"
                       )

    parser.add_argument( "--comoving", type=str2bool, metavar="BOOLEAN", gamer_name="COMOVING",
                         default=False,
                         depend={"gravity":True},
                         constraint={ True:{"eos":"GAMMA"} },
                         help="Comoving frame for cosmological simulations.\n"
                       )

    # A.4 particle
    parser.add_argument( "--particle", type=str2bool, metavar="BOOLEAN", gamer_name="PARTICLE",
                         default=False,
                         help="Enable particles.\n"
                       )
    parser.add_argument( "--tracer", type=str2bool, metavar="BOOLEAN", gamer_name="TRACER",
                         default=False,
                         depend={"particle":True},
                         help="Enable tracer particles.\n"
                       )

    parser.add_argument( "--store_par_acc", type=str2bool, metavar="BOOLEAN", gamer_name="STORE_PAR_ACC",
                         default=True,
                         depend={"particle":True},
                         help="Store particle acceleration (recommended).\n"
                       )

    parser.add_argument( "--star_formation", type=str2bool, metavar="BOOLEAN", gamer_name="STAR_FORMATION",
                         default=False,
                         depend={"particle":True},
                         help="Allow creating new particles after initialization. Must trun on <--store_pot_ghost> when <--store_par_acc> is adoped.\n"
                       )

    parser.add_argument( "--feedback", type=str2bool, metavar="BOOLEAN", gamer_name="FEEDBACK",
                         default=False,
                         depend={"particle":True},
                         help="Feedback from particles to grids and vice versa.\n"
                       )

    parser.add_argument( "--par_attribute", type=int, metavar="INTEGER", gamer_name="PAR_NATT_USER",
                         default=0,
                         depend={"particle":True},
                         help="Set the number of user-defined particle attributes.\n"
                       )

    # A.5 grackle
    parser.add_argument( "--grackle", type=str2bool, metavar="BOOLEAN", gamer_name="SUPPORT_GRACKLE",
                         default=False,
                         constraint={ True:{"model":"HYDRO", "eos":"GAMMA"} },
                         help="Enable Grackle, a chemistry and radiative cooling library. Must set <--passive> according to the primordial chemistry network set by GRACKLE_PRIMORDIAL. Please enable OpenMP when compiling Grackle (by 'make omp-on').\n"
                       )

    # B. miscellaneous options
    parser.add_argument( "--nlevel", type=int, metavar="INTEGER", gamer_name="NLEVEL",
                         default=10,
                         help="Set the total number of AMR levels including the root level.\n"
                       )

    parser.add_argument( "--max_patch", type=int, metavar="INTEGER", gamer_name="MAX_PATCH",
                         default=1000000,
                         help="Set the maximum number of patches on each AMR level.\n"
                       )

    parser.add_argument( "--patch_size", type=int, metavar="INTEGER", gamer_name="PATCH_SIZE",
                         default=8,
                         help="Set the number of cells along each direction in a single patch. Must be an even number greater than or equal to 8.\n"
                       )

    parser.add_argument( "--debug", type=str2bool, metavar="BOOLEAN", gamer_name="GAMER_DEBUG",
                         default=False,
                         help="Enable debug mode.\n"
                       )

    parser.add_argument( "--bitwise_reproducibility", type=str2bool, metavar="BOOLEAN", gamer_name="BITWISE_REPRODUCIBILITY",
                         default=None,
                         help="Enable bitwise reproducibility.\n"
                       )

    parser.add_argument( "--timing", type=str2bool, metavar="BOOLEAN", gamer_name="TIMING",
                         default=True,
                         help="Enable timing analysis of a simulation.\n"
                       )

    parser.add_argument( "--timing_solver", type=str2bool, metavar="BOOLEAN", gamer_name="TIMING_SOLVER",
                         default=False,
                         constraint={ True:{"timing":True} },
                         help="Enable more detailed timing analysis of GPU solvers. This option will disable CPU/GPU overlapping and thus deteriorate performance. Must enable <--timing>.\n"
                       )

    parser.add_argument( "--double", type=str2bool, metavar="BOOLEAN", gamer_name="FLOAT8",
                         default=False,
                         help="Enable double precision.\n"
                       )

    parser.add_argument( "--laohu", type=str2bool, metavar="BOOLEAN", gamer_name="LAOHU",
                         default=False,
                         help="Work on the NAOC Laohu GPU cluster.\n"
                       )

    parser.add_argument( "--hdf5", type=str2bool, metavar="BOOLEAN", gamer_name="SUPPORT_HDF5",
                         default=False,
                         help="Support HDF5 (recommended).\n"
                       )

    parser.add_argument( "--gsl", type=str2bool, metavar="BOOLEAN", gamer_name="SUPPORT_GSL",
                         default=False,
                         help="Support GNU scientific library.\n"
                       )

    parser.add_argument( "--fftw", type=str, metavar="TYPE", gamer_name="SUPPORT_FFTW",
                         default=NONE_STR, choices=[NONE_STR, "FFTW2", "FFTW3"],
                         help="Support FFTW library.\n"
                       )

    parser.add_argument( "--libyt", type=str2bool, metavar="BOOLEAN", gamer_name="SUPPORT_LIBYT",
                         default=False,
                         help="Support yt inline analysis.\n"
                       )

    parser.add_argument( "--libyt_patchgroup", type=str2bool, metavar="BOOLEAN", gamer_name="LIBYT_USE_PATCH_GROUP",
                         default=True,
                         depend={"libyt":True},
                         help="Use patch groups instead of patches as the unit in libyt for better performance (recommended). Must enable <--libyt>.\n"
                       )

    parser.add_argument( "--libyt_interactive", type=str2bool, metavar="BOOLEAN", gamer_name="LIBYT_INTERACTIVE",
                         default=False,
                         depend={"libyt":True},
                         help="Enable the interactive mode of libyt. This activates python prompt and does not shut down a simulation when there are errors in an inline python script. Must compile libyt with INTERACTIVE_MODE. Must enable <--libyt>.\n"
                       )

    parser.add_argument( "--rng", type=str, metavar="TYPE", gamer_name="RANDOM_NUMBER",
                         default="RNG_GNU_EXT",
                         choices=["RNG_GNU_EXT", "RNG_CPP11"],
                         help="Select the random number generator (RNG_GNU_EXT: GNU extension drand48_r, RNG_CPP11: c++11 <random>).\nRNG_GNU_EXT may not be supported on some macOS.\nFor RNG_CPP11, add -std=c++11 to CXXFLAG in your config file.\n"
                       )

    # C. parallelization and flags
    parser.add_argument( "--openmp", type=str2bool, metavar="BOOLEAN", gamer_name="OPENMP",
                         default=True,
                         help="Enable OpenMP parallization.\n"
                       )

    parser.add_argument( "--mpi", type=str2bool, metavar="BOOLEAN", gamer_name="LOAD_BALANCE=HILBERT",
                         default=False,
                         help="Enable MPI parallelization.\n"
                       )

    parser.add_argument( "--overlap_mpi", type=str2bool, metavar="BOOLEAN", gamer_name="OVERLAP_MPI",
                         default=False,
                         constraint={ True:{"mpi":True} },
                         help="Overlap MPI communication with computation. Not supported yet. Must enable <--mpi>.\n"
                       )

    parser.add_argument( "--gpu", type=str2bool, metavar="BOOLEAN", gamer_name="GPU",
                         default=False,
                         help="Enable GPU. Must set <--gpu_arch> as well.\n"
                       )

    parser.add_argument( "--gpu_arch", type=str, metavar="TYPE", gamer_name="GPU_ARCH",
                         depend={"gpu":True},
                         default="TURING", choices=["FERMI", "KEPLER", "MAXWELL", "PASCAL", "VOLTA", "TURING", "AMPERE"],
                         help="Select the architecture of GPU.\n"
                       )

    args, name_table, depends, constraints = parser.parse_args()
    args = vars( args )

    # 1. Print out a detailed help message then exit.
    if args["lh"]:
        parser.print_help_detail()
        exit()

    # 2. Conditional default arguments.
    args = set_conditional_defaults( args )
    return args, name_table, depends, constraints

def load_config( config ):
    print("Using %s as the config."%(config))
    paths, compilers, flags = {}, {"CXX":"", "CXX_MPI":""}, {"CXXFLAG":"", "OPENMPFLAG":"", "LIBFLAG":"", "CUDAFLAG":""}
    with open( config, 'r') as f:
        lines = f.readlines()

    for line in lines:
        temp = list( filter( None, re.split(" |:=|\n", line) ) ) # separate by " " and ":="
        if len(temp) == 0: continue             # empty line
        if temp[0][0] == "#": continue          # skip comment line
        if temp[0] in flags:
            if len(temp) == 1: continue         # empty flag
            for i in range(1, len(temp)):
                if temp[i][0] == "#": break     # commented out
                flags[temp[0]] += temp[i] + " "
        elif temp[0] in compilers:
            if len(temp) == 1: continue         # empty compiler
            if temp[1][0] == "#": continue      # comment out
            if compilers[temp[0]] != "": color_print("Warning: The original compiler will be overwritten. <%s>: %s --> %s"%(temp[0], compilers[temp[0]], temp[1]), BCOLOR.WARNING)
            compilers[temp[0]] = temp[1]
        else:
            try:
               paths[temp[0]] = temp[1]
            except:
               paths[temp[0]] = ''

    return paths, compilers, flags

def set_conditional_defaults( args ):
    if args["unsplit_gravity"] == None:
        args["unsplit_gravity"] = True if args["model"] == "HYDRO" else False

    if args["bitwise_reproducibility"] == None:
        args["bitwise_reproducibility"] = True if args["debug"] else False

    return args

def set_sims( name_table, depends, **kwargs ):
    opt_str = ""
    # Loop all the simulation options in GAMER.
    for opt, gamer_name in name_table.items():
        store = True
        # check if the depend is true
        if opt in depends:
            for depend, val in depends[opt].items():
                if type(val) != type([]): val = [val]   # transform to list
                if kwargs[depend] not in val: store = False

        if not store: continue
        if opt == "eos":        # special string prefix of EOS
            opt_str = add_option( opt_str, name=gamer_name, val="EOS_"+kwargs[opt] )
        else:
            opt_str = add_option( opt_str, name=gamer_name, val=kwargs[opt] )

    # Hard-code the option of serial.
    if not kwargs["mpi"]: opt_str = add_option( opt_str, name="SERIAL", val=True )

    return {"SIMU_OPTION":opt_str}

def set_compile( paths, compilers, flags, kwargs ):
    com_opt = {}

    # 1. Set the complier.
    com_opt["CXX"] = paths["MPI_PATH"]+"/bin/"+compilers["CXX_MPI"] if kwargs["mpi"] else compilers["CXX"]

    # 2. Set the OpenMP flags.
    if not kwargs["openmp"]: flags["OPENMPFLAG"] = ""

    # 3. Write flags to complie option dictionary.
    for key, val in flags.items():
        com_opt[key] = val

    return com_opt

def validation( paths, depends, constraints, **kwargs ):
    success = True

    # 0. Checking the Makefile_base existance.
    if not os.path.isfile( GAMER_MAKE_BASE ):
        color_print("ERROR: %s does not exist."%(GAMER_MAKE_BASE), BCOLOR.FAIL)
        success = False

    # 1. Checking general constraints.
    for opt, condition in constraints.items():
        check = True
        if opt in depends:
            for depend, val in depends[opt].items():
                if type(val) != type([]): val = [val]   # transform to list
                if kwargs[depend] not in val: check = False

        if not check: continue  # do not validate if the depend is not set

        for opt_val, check in condition.items():
            if kwargs[opt] != opt_val: continue         # not the value to be checked
            for check_opt, check_val in check.items():
                if type(check_val) != type([]): check_val = [check_val]   # transform to list
                if kwargs[check_opt] in check_val: continue     # satisify the validation

                val_str = ', '.join(str(x) for x in check_val)
                color_print("ERROR: The option <--%s=%s> requires <--%s> to be set to [%s]. Current: <--%s=%s>."%(opt, str(kwargs[opt]), check_opt, val_str, check_opt, kwargs[check_opt]), BCOLOR.FAIL)
                success = False

    # 2. Checking other conditions.
    # A. Physics
    # A.1 Module
    if kwargs["model"] == "HYDRO":
        if kwargs["passive"] < 0:
            color_print("ERROR: Passive scalar should not be negative. Current: %d"%kwargs["passive"], BCOLOR.FAIL)
            success = False

        if kwargs["dual"] not in [NONE_STR, "DE_ENPY"]:
            color_print("ERROR: This dual energy form is not supported yet. Current: %s"%kwargs["dual"], BCOLOR.FAIL)
            success = False

        if kwargs["cosmic_ray"]:
            color_print("ERROR: <--cosmic_ray> is not supported yet.", BCOLOR.FAIL)
            success = False

    elif kwargs["model"] == "ELBDM":
        if kwargs["passive"] < 0:
            color_print("ERROR: Passive scalar should not be negative. Current: %d"%kwargs["passive"], BCOLOR.FAIL)
            success = False

    elif kwargs["model"] == "PAR_ONLY":
        color_print("ERROR: <--model=PAR_ONLY> is not supported yet.", BCOLOR.FAIL)
        success = False
    else:
        color_print("ERROR: Unrecognized model: <%s>. Please add to the model choices."%kwargs["model"], BCOLOR.FAIL)
        success = False

    # A.3 Particle
    if kwargs["particle"]:
        if kwargs["star_formation"] and kwargs["store_par_acc"] and not kwargs["store_pot_ghost"]:
            color_print("ERROR: <--store_pot_ghost> must be enabled when <--star_formation> and <--store_par_acc> are enabled.", BCOLOR.FAIL)
            success = False
        if not kwargs["gravity"] and not kwargs["tracer"]:
            color_print("ERROR: At least one of <--gravity> or <--tracer> must be enabled for <--particle>.", BCOLOR.FAIL)
            success = False
        if kwargs["par_attribute"] < 0:
            color_print("ERROR: Number of particle attributes should not be negative. Current: %d"%kwargs["par_attribute"], BCOLOR.FAIL)
            success = False

    # B. miscellaneous options
    if kwargs["nlevel"] < 1:
        color_print("ERROR: <--nlevel> should be greater than zero. Current: %d"%kwargs["nlevel"], BCOLOR.FAIL)
        success = False

    if kwargs["max_patch"] < 1:
        color_print("ERROR: <--max_patch> should be greater than zero. Current: %d"%kwargs["max_patch"], BCOLOR.FAIL)
        success = False

    if kwargs["patch_size"]%2 != 0 or kwargs["patch_size"] < 8:
        color_print("ERROR: <--patch_size> should be an even number greater than or equal to 8. Current: %d"%kwargs["patch_size"], BCOLOR.FAIL)
        success = False

    if kwargs["overlap_mpi"]:
        color_print("ERROR: <--overlap_mpi> is not supported yet.", BCOLOR.FAIL)
        success = False

    if not success: raise BaseException(BCOLOR.FAIL+"The above validation failed."+BCOLOR.ENDC)
    return


def warning( paths, **kwargs ):
    # 1. Makefile
    if os.path.isfile( GAMER_MAKE_OUT ):
        color_print("Warning: %s already exists and will be overwritten."%(GAMER_MAKE_OUT), BCOLOR.WARNING)

    # 2. Physics
    if kwargs["model"] == "ELBDM" and kwargs["passive"] != 0:
        color_print("Warning: Not supported yet and can only be used as auxiliary fields.", BCOLOR.WARNING)

    # 3. Path
    if kwargs["gpu"]:
        if paths.setdefault("CUDA_PATH", "") == "":
            color_print("Warning: CUDA_PATH is not given in %s.config when enabling <--gpu>."%(kwargs["machine"]), BCOLOR.WARNING)

    if kwargs["fftw"] == "FFTW2":
        if paths.setdefault("FFTW2_PATH", "") == "":
            color_print("Warning: FFTW2_PATH is not given in %s.config when enabling <--fftw=FFTW2>."%(kwargs["machine"]), BCOLOR.WARNING)

    if kwargs["fftw"] == "FFTW3":
        if paths.setdefault("FFTW3_PATH", "") == "":
            color_print("Warning: FFTW3_PATH is not given in %s.config when enabling <--fftw=FFTW3>."%(kwargs["machine"]), BCOLOR.WARNING)

    if kwargs["mpi"]:
        if paths.setdefault("MPI_PATH", "") == "":
            color_print("Warning: MPI_PATH is not given in %s.config when enabling <--mpi>."%(kwargs["machine"]), BCOLOR.WARNING)

    if kwargs["hdf5"]:
        if paths.setdefault("HDF5_PATH", "") == "":
            color_print("Warning: HDF5_PATH is not given in %s.config when enabling <--hdf5>."%(kwargs["machine"]), BCOLOR.WARNING)

    if kwargs["grackle"]:
        if paths.setdefault("GRACKLE_PATH", "") == "":
            color_print("Warning: GRACKLE_PATH is not given in %s.config when enabling <--grackle>."%(kwargs["machine"]), BCOLOR.WARNING)

    if kwargs["gsl"]:
        if paths.setdefault("GSL_PATH", "") == "":
            color_print("Warning: GSL_PATH is not given in %s.config when enabling <--gsl>."%(kwargs["machine"]), BCOLOR.WARNING)

    if kwargs["libyt"]:
        if paths.setdefault("LIBYT_PATH", "") == "":
            color_print("Warning: LIBYT_PATH is not given in %s.config when enabling <--libyt>."%(kwargs["machine"]), BCOLOR.WARNING)

    return



####################################################################################################
# Main execution
####################################################################################################
# 1. Load the input arguments
args, name_table, depends, constraints = load_arguments()

#------------------------------------------------------------
# 2. Prepare the makefile args
# 2.1 Load the machine setup
paths, compilers, flags = load_config( "%s/%s.config"%(GAMER_CONFIG_DIR, args["machine"]) )

# 2.2 Validate arguments
validation( paths, depends, constraints, **args )

warning( paths, **args )

# 2.3 add the SIMU_OPTION
print("")
print("========================================")
print("GAMER has the following setting.")
print("----------------------------------------")
sims = set_sims( name_table, depends, **args )

# 2.4 setup compiler
compiles = set_compile( paths, compilers, flags, args )


#------------------------------------------------------------
# 3. Create Makefile
# 3.1 Read
with open( GAMER_MAKE_BASE, "r" ) as make_base:
    makefile = make_base.read()

# 3.2 Replace
print("----------------------------------------")
for key, val in paths.items():
    print("%-15s : %s"%(key, val))
    makefile, num = re.subn(r"@@@%s@@@"%(key), val, makefile)
    if num == 0: raise BaseException("The string @@@%s@@@ is not replaced correctly."%key)

for key, val in sims.items():
    makefile, num = re.subn(r"@@@%s@@@"%(key), val, makefile)
    if num == 0: raise BaseException("The string @@@%s@@@ is not replaced correctly."%key)

print("----------------------------------------")
for key, val in compiles.items():
    print("%-10s : %s"%(key, val))
    makefile, num = re.subn(r"@@@%s@@@"%(key), val, makefile)
    if num == 0: raise BaseException("The string @@@%s@@@ is not replaced correctly."%key)

# 3.3 Write
with open( GAMER_MAKE_OUT, "w") as make_out:
    make_out.write( makefile )


print("========================================")
print("%s is created."%GAMER_MAKE_OUT)
print("========================================")

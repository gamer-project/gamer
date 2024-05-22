#!/usr/bin/python3
"""
User and developer guides of this script are provided in the following link.

   https://github.com/gamer-project/gamer/wiki/Installation%3A-Configure.py

"""

####################################################################################################
# Packages
####################################################################################################
import argparse
import logging
import os
import sys
import re
import ctypes



####################################################################################################
# Global variables
####################################################################################################
NONE_STR = "OFF"

CLOSE_DIST  = 2
PRINT_WIDTH = 100

GAMER_CONFIG_DIR  = os.path.join("..", "configs")
GAMER_MAKE_BASE   = "Makefile_base"
GAMER_MAKE_OUT    = "Makefile"
GAMER_DESCRIPTION = "Prepare a customized Makefile for GAMER.\nDefault values are marked by '*'.\nUse -lh to show a detailed help message.\n"
GAMER_EPILOG      = "2023 Computational Astrophysics Lab, NTU. All rights reserved.\n"

LOGGER     = logging.getLogger()
LOG_FORMAT = '%(asctime)s %(levelname)-8s: %(message)s'
logging.basicConfig( filename=GAMER_MAKE_OUT+'.log', filemode='w', level=logging.INFO, format=LOG_FORMAT )



####################################################################################################
# Classes
####################################################################################################
class CustomFormatter( logging.Formatter ):
    """
    See: https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
    """
    HEADER    = '\033[95m'
    OKBLUE    = '\033[94m'
    OKCYAN    = '\033[96m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'

    CONSOLE_FORMAT = '%(levelname)-8s: %(message)s'
    FORMATS = { logging.DEBUG   : HEADER  + CONSOLE_FORMAT + ENDC,
                logging.INFO    : ENDC    + "%(message)s"  + ENDC,
                logging.WARNING : WARNING + CONSOLE_FORMAT + ENDC,
                logging.ERROR   : FAIL    + CONSOLE_FORMAT + ENDC,
                logging.CRITICAL: FAIL    + CONSOLE_FORMAT + ENDC }

    def format( self, record ):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

class ArgumentParser( argparse.ArgumentParser ):
    def __init__( self, *args, **kwargs ):
        self.program     = kwargs.copy()
        self.options     = []
        self.depends     = {}
        self.constraints = {}
        self.gamer_names = {}
        super(ArgumentParser, self).__init__(*args, **kwargs)

    def add_argument( self, *args, **kwargs ):
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
        option = kwargs.copy()
        option["flags"] = [ item for item in args ]
        self.options.append(option)

    def parse_args( self, args=None, namespace=None ):
        args, argv = self.parse_known_args(args, namespace)
        msg = "\n"
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
            if min_dist <= CLOSE_DIST: msg += ', do you mean: %s ?\n'%(pos_key)
            msg += '\n'
            if arg == '--gpu_arch': msg += "ERROR: <--gpu_arch> is deprecated. Please set <GPU_COMPUTE_CAPABILITY> in your machine *.config file (see ../configs/template.config).\n"

        if len(argv) != 0: self.error( msg )
        return args, self.gamer_names, self.depends, self.constraints

    def _get_option_tuples(self, option_string):
        # This function is directly from the source code of `argparse`.
        # We decided to add the function manually because versions prior to Python 3.5 do not support `allow_abbrev`.
        # See: https://github.com/python/cpython/blob/main/Lib/argparse.py
        result = []

        # option strings starting with two prefix characters are only split at the '='
        chars = self.prefix_chars
        if option_string[0] in chars and option_string[1] in chars:
            pass # we always use `allow_abbrev=False`

        # single character options can be concatenated with their arguments
        # but multiple character options always have to have their arguments separate
        elif option_string[0] in chars and option_string[1] not in chars:
            option_prefix = option_string
            short_option_prefix = option_string[:2]
            short_explicit_arg = option_string[2:]

            for option_string in self._option_string_actions:
                if option_string == short_option_prefix:
                    action = self._option_string_actions[option_string]
                    tup = action, option_string, '', short_explicit_arg
                    result.append(tup)
                elif option_string.startswith(option_prefix):
                    action = self._option_string_actions[option_string]
                    tup = action, option_string, None, None
                    result.append(tup)

        # shouldn't ever get here
        else:
            self.error(_('unexpected option string: %s') % option_string)

        return result # return the collected option tuples

    def print_usage( self, *args, **kwargs ):
        if "usage" in self.program:
            print("Usage: %s\n" % self.program["usage"])
            return

        usage = []
        for option in self.options:
            for item in option["flags"]:
                # the order of if does matter here
                possibles = ""
                if   "choices" in option: possibles += "{%s}"%(", ".join([ str(opt) for opt in option["choices"] ]))
                elif "metavar" in option: possibles += option["metavar"]
                elif "dest"    in option: possibles += option["dest"] if "default" in option else option["dest"].upper()
                else:                     possibles += re.sub(r"^(-{1,})", "", item).upper()

                default_value = ""
                if "default" in option:
                    default_value += "*"
                    default_value += "Depend" if option["default"] is None else str(option["default"])

                if "action" in option:
                    if option["action"] in ["help", "store_const", "store_true", "store_false"]:
                        possibles = "\b"
                        default_value = "\b"

                usage += [ "[%s]"%(" ".join([item, possibles, default_value])) ]

        indent = "Usage: %s " % os.path.basename(sys.argv[0])
        output = indent + " " + str.join(" ", usage) + "\n"
        print( string_align(output, indent, PRINT_WIDTH, "]") )

    def print_option( self ):
        print("Options:")
        option_indent = 0
        for option in self.options:
            option["flags2"] = str.join(", ", [ "%s %s" % (item, option["metavar"]) if "metavar" in option else "%s %s" % (item, option["dest"].upper()) if "dest" in option else item for item in option["flags"] ])
            if len(option["flags2"]) > option_indent:
                option_indent = len(option["flags2"])

        template = "  %-" + str(option_indent) + "s  "
        for option in self.options:
            indent = template %(option["flags2"])
            output = indent

            if "help" in option: output += option["help"]

            if "action" in option:
                if option["action"] == "help":
                    print( string_align(output, indent, PRINT_WIDTH, " ") )
                    continue

            if "choices" in option:
                temp = [ str(opt) for opt in option["choices"] ]
                output += "Choice: [%s] => "%(", ".join(temp))

            if "default" in option:
                output += "Default: %s" %("Depend" if option["default"] is None else str(option["default"]))

            if "action" in option:
                output += "Default: False" if option["action"] == "store_true" else "Default: True"

            print( string_align(output, indent, PRINT_WIDTH, " ") )

    def print_help( self, *args, **kwargs ):
        # print usage, description, options, then epilog
        self.print_usage()
        if "description"  in self.program: print(self.program["description"])
        if "print_detail" in kwargs: self.print_option()
        if "epilog"       in self.program: print(self.program["epilog"])



####################################################################################################
# Functions
####################################################################################################
def str2bool( v ):
    if isinstance(v, bool): return v

    if   v.lower() == "true":  return True
    elif v.lower() == "false": return False
    else: raise TypeError("Can not convert <%s> to boolean."%(v))
    return

def add_option( opt_str, name, val ):
    # NOTE: Every -Doption must have a trailing space.
    if type(val) == type(True):
        if val: opt_str += "-D%s "%(name)
        LOGGER.info("%-25s : %r"%(name, val))
    elif type(val) == type("str"):
        if val != NONE_STR:
            opt_str += "-D%s=%s "%(name, val)
            LOGGER.info("%-25s : %s"%(name, val))
    elif type(val) == type(0):
        opt_str += "-D%s=%d "%(name, val)
        LOGGER.info("%-25s : %d"%(name, val))
    elif type(val) == type(0.):
        opt_str += "-D%s=%f "%(name, val)
        LOGGER.info("%-25s : %f"%(name, val))
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

def get_gpu_compute_capability():
    """
    Outputs some information on CUDA-enabled devices on your computer, including current memory usage.

    It's a port of https://gist.github.com/f0k/0d6431e3faa60bffc788f8b4daa029b1
    from C to Python with ctypes, so it can run without compiling anything. Note
    that this is a direct translation with no attempt to make the code Pythonic.
    It's meant as a general demonstration on how to obtain CUDA device information
    from Python without resorting to nvidia-smi or a compiled Python extension.

    Author: Jan Schluter
    License: MIT (https://gist.github.com/f0k/63a664160d016a491b2cbea15913d549#gistcomment-3870498)
    Others: https://en.wikipedia.org/wiki/CUDA#GPUs_supported
    """
    CUDA_SUCCESS = 0
    libnames = ('libcuda.so', 'libcuda.dylib', 'cuda.dll')
    for libname in libnames:
        try:
            cuda = ctypes.CDLL(libname)
        except OSError:
            continue
        else:
            break
    else:
        raise OSError("could not load any of: " + ' '.join(libnames))

    nGpus, cc_major, cc_minor, device = ctypes.c_int(), ctypes.c_int(), ctypes.c_int(), ctypes.c_int()

    def cuda_check_error( result ):
        if result == CUDA_SUCCESS: return

        error_str = ctypes.c_char_p()

        cuda.cuGetErrorString(result, ctypes.byref(error_str))
        raise BaseException( "CUDA failed with error code %d: %s"%( result, error_str.value.decode() ) )

        return

    cuda_check_error( cuda.cuInit(0) )
    cuda_check_error( cuda.cuDeviceGetCount(ctypes.byref(nGpus)) )

    if nGpus.value > 1: LOGGER.warning("More than one GPU --> select the compute capability of the last GPU.")
    for i in range(nGpus.value):
        cuda_check_error( cuda.cuDeviceGet(ctypes.byref(device), i) )
        cuda_check_error( cuda.cuDeviceComputeCapability(ctypes.byref(cc_major), ctypes.byref(cc_minor), device) )

    compute_capability = cc_major.value*100 + cc_minor.value*10
    return compute_capability

def string_align( string, indent_str, width, end_char ):
    """
    end_char : The ending character of a word.
    """
    N          = len(indent_str)
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

def load_arguments():
    parser = ArgumentParser( description = GAMER_DESCRIPTION,
                             formatter_class = argparse.RawTextHelpFormatter,
                             epilog = GAMER_EPILOG,
                             add_help = False
                           )

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

    # A. options of diffierent physical models
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
                         help="The hydrodynamic/MHD integrator. MHD only supports MHM, MHM_RP and CTU.\n"
                       )

    parser.add_argument( "--slope", type=str, metavar="TYPE", gamer_name="LR_SCHEME",
                         default="PPM", choices=["PLM", "PPM"],
                         depend={"model":"HYDRO"},
                         help="The spatial data reconstruction method (PLM: piecewise-linear, PPM: piecewise-parabolic). Useless for <--flu_scheme=RTVD>.\n"
                       )

    parser.add_argument( "--flux", type=str, metavar="TYPE", gamer_name="RSOLVER",
                         default=None, choices=["EXACT", "ROE", "HLLE", "HLLC", "HLLD"],
                         depend={"model":"HYDRO"},
                         constraint={ "ROE":{"eos":"GAMMA"},
                                      "EXACT":{"eos":"GAMMA"} },
                         help="The Riemann solver. Pure hydro: EXACT/ROE/HLLE/HLLC^, MHD: ROE/HLLE/HLLD^, SRHD: HLLE/HLLC^, (^ indicates the recommended and default solvers). Useless for RTVD.\n"
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
                         constraint={ True:{"flu_scheme":["MHM", "MHM_RP", "CTU"], "flux":["ROE", "HLLE", "HLLD"]},
                                     False:{"flux":["EXACT", "ROE", "HLLE", "HLLC"]} },
                         help="Magnetohydrodynamics.\n"
                       )

    parser.add_argument( "--srhd", type=str2bool, metavar="BOOLEAN", gamer_name="SRHD",
                         default=False,
                         depend={"model":"HYDRO"},
                         constraint={ True:{"flu_scheme":["MHM", "MHM_RP"], "flux":["HLLE", "HLLC"], "eos":["TAUBMATHEWS"], "dual":[NONE_STR], "mhd":False, "gravity":False} },
                         help="Special Relativistic Hydrodynamics.\n"
                       )

    parser.add_argument( "--cosmic_ray", type=str2bool, metavar="BOOLEAN", gamer_name="COSMIC_RAY",
                         default=False,
                         depend={"model":"HYDRO"},
                         constraint={ True:{"dual":[NONE_STR], "eos":"COSMIC_RAY", "comoving":False} },
                         help="Enable cosmic ray. Must use <--eos=COSMIC_RAY>.\n"
                       )

    parser.add_argument( "--eos", type=str, metavar="TYPE", gamer_name="EOS",
                         default=None, choices=["GAMMA", "ISOTHERMAL", "NUCLEAR", "TABULAR", "COSMIC_RAY", "TAUBMATHEWS", "USER"],
                         depend={"model":"HYDRO"},
                         constraint={ "ISOTHERMAL":{"barotropic":True}, "COSMIC_RAY":{"cosmic_ray":True}, "TAUBMATHEWS":{"srhd":True} },
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
                         help="Use unsplitting method to couple gravity to the target model (recommended). Supported only for <--model=HYDRO>. When <--model=HYDRO>, the default is True.\n"
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

    parser.add_argument( "--double_par", type=str2bool, metavar="BOOLEAN", gamer_name="FLOAT8_PAR",
                         default=None,
                         depend={"particle":True},
                         help="Enable double precision for particle attributes.\n"
                       )
    # A.5 grackle
    parser.add_argument( "--grackle", type=str2bool, metavar="BOOLEAN", gamer_name="SUPPORT_GRACKLE",
                         default=False,
                         constraint={ True:{"model":"HYDRO", "eos":["GAMMA", "COSMIC_RAY"], "comoving":False} },
                         help="Enable Grackle, a chemistry and radiative cooling library. Must set <--passive> according to the primordial chemistry network set by GRACKLE_PRIMORDIAL. Please enable OpenMP when compiling Grackle (by 'make omp-on').\n"
                       )

    # A.6 microphysics
    parser.add_argument( "--cr_diffusion", type=str2bool, metavar="BOOLEAN", gamer_name="CR_DIFFUSION",
                         default=False,
                         depend={"cosmic_ray":True},
                         constraint={ True:{"cosmic_ray":True, "mhd":True} },
                         help="Enable cosmic-ray diffusion. Must enable <--mhd> and <--cosmic_ray>.\n"
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
                         help="Enable bitwise reproducibility. When <--debug=True>, the default is True.\n"
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

    parser.add_argument( "--libyt_reload", type=str2bool, metavar="BOOLEAN", gamer_name="LIBYT_RELOAD",
                         default=False,
                         depend={"libyt":True},
                         help="Allow for reloading libyt scripts during runtime. Must compile libyt with INTERACTIVE_MODE. Must enable <--libyt>.\n"
                       )

    parser.add_argument( "--libyt_jupyter", type=str2bool, metavar="BOOLEAN", gamer_name="LIBYT_JUPYTER",
                         default=False,
                         depend={"libyt":True},
                         help="Allow for in situ analysis using Jupyter Notebook / JupyterLab through libyt. Must compile libyt with JUPYTER_KERNEL. Must enable <--libyt>.\n"
                       )

    parser.add_argument( "--rng", type=str, metavar="TYPE", gamer_name="RANDOM_NUMBER",
                         default="RNG_GNU_EXT",
                         choices=["RNG_GNU_EXT", "RNG_CPP11"],
                         help="Select the random number generator (RNG_GNU_EXT: GNU extension drand48_r, RNG_CPP11: c++11 <random>).\nRNG_GNU_EXT may not be supported on some macOS.\nFor RNG_CPP11, add -std=c++11 to CXXFLAG in your config file.\n"
                       )

    # C. parallelization and flags
    parser.add_argument( "--openmp", type=str2bool, metavar="BOOLEAN", gamer_name="OPENMP",
                         default=True,
                         help="Enable OpenMP parallelization.\n"
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
                         help="Enable GPU. Must set <GPU_COMPUTE_CAPABILITY> in your machine *.config file as well.\n"
                       )

    args, name_table, depends, constraints = parser.parse_args()
    args = vars( args )

    # 1. Print out a detailed help message then exit.
    if args["lh"]:
        parser.print_help(print_detail=True)
        exit()

    # 2. Conditional default arguments.
    args = set_conditional_defaults( args )
    return args, name_table, depends, constraints

def load_config( config ):
    LOGGER.info("Using %s as the config."%(config))
    paths, compilers = {}, {"CXX":"", "CXX_MPI":""}
    flags = {"CXXFLAG":"", "OPENMPFLAG":"", "LIBFLAG":"", "NVCCFLAG_COM":"", "NVCCFLAG_FLU":"", "NVCCFLAG_POT":""}
    gpus  = {"GPU_COMPUTE_CAPABILITY":""}
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
            if temp[1][0] == "#": continue      # commented out
            if compilers[temp[0]] != "": LOGGER.warning("The original compiler will be overwritten. <%s>: %s --> %s"%(temp[0], compilers[temp[0]], temp[1]))
            compilers[temp[0]] = temp[1]
        elif temp[0] in gpus:
            if len(temp) == 1: continue         # empty
            if temp[1][0] == "#": continue      # commented out
            if gpus[temp[0]] != "": LOGGER.warning("The original value will be overwritten. <%s>: %s --> %s"%(temp[0], gpus[temp[0]], temp[1]))
            gpus[temp[0]] = temp[1]
        else:
            try:
               paths[temp[0]] = temp[1]
            except:
               paths[temp[0]] = ''

    return paths, compilers, flags, gpus

def set_conditional_defaults( args ):
    if args["unsplit_gravity"] is None:
        args["unsplit_gravity"] = (args["model"] == "HYDRO")

    if args["bitwise_reproducibility"] is None:
        args["bitwise_reproducibility"] = args["debug"]

    if args["double_par"] is None:
        args["double_par"] = args["double"]

    if args["flux"] is None:
        args["flux"] = "HLLD" if args["mhd"] else "HLLC"

    if args["eos"] is None:
        if   args["cosmic_ray"]: args["eos"] = "COSMIC_RAY"
        elif args["srhd"]      : args["eos"] = "TAUBMATHEWS"
        else                   : args["eos"] = "GAMMA"
    return args

def set_gpu( gpus, flags, args ):
    gpu_opts = {}
    compute_capability = gpus["GPU_COMPUTE_CAPABILITY"]

    # 1. Check the compute capability
    if compute_capability == "":
        if args["gpu"]: raise ValueError("GPU_COMPUTE_CAPABILITY is not set in `../configs/%s.config`. See `../configs/template.config` for illustration."%args["machine"])
        return gpu_opts
    compute_capability = int(compute_capability)

    if   compute_capability < 0:
        compute_capability = get_gpu_compute_capability()
    elif compute_capability < 200:
        raise ValueError("Incorrect GPU_COMPUTE_CAPABILITY range (>=200)")
    gpu_opts["GPU_COMPUTE_CAPABILITY"] = str(compute_capability)

    # 2. Set NVCCFLAG_ARCH
    flag_num = compute_capability // 10
    gpu_opts["NVCCFLAG_ARCH"] = '-gencode arch=compute_%d,code=\\"compute_%d,sm_%d\\"'%(flag_num, flag_num, flag_num)

    # 3. Set MAXRREGCOUNT_FLU
    if 300 <= compute_capability and compute_capability <= 370:
        if args["double"]:
            gpu_opts["MAXRREGCOUNT_FLU"] = "--maxrregcount=128"
        else:
            gpu_opts["MAXRREGCOUNT_FLU"] = "--maxrregcount=70"
    elif 500 <= compute_capability and compute_capability <= 870:
        if args["double"]:
            gpu_opts["MAXRREGCOUNT_FLU"] = "--maxrregcount=192"
        else:
            gpu_opts["MAXRREGCOUNT_FLU"] = "--maxrregcount=128"
    return gpu_opts

def set_sims( name_table, depends, **kwargs ):
    opt_str = ""
    # loop all the simulation options in GAMER.
    for opt, gamer_name in name_table.items():
        store = True
        # check if depend is true
        if opt in depends:
            for depend, val in depends[opt].items():
                if type(val) != type([]): val = [val]   # transform to list
                if kwargs[depend] not in val: store = False

        if not store: continue
        if opt == "eos":        # special string prefix of EOS
            opt_str = add_option( opt_str, name=gamer_name, val="EOS_"+kwargs[opt] )
        else:
            opt_str = add_option( opt_str, name=gamer_name, val=kwargs[opt] )

    # hard-code the option of serial.
    if not kwargs["mpi"]: opt_str = add_option( opt_str, name="SERIAL", val=True )

    return {"SIMU_OPTION":opt_str}

def set_compile( paths, compilers, flags, kwargs ):
    com_opt = {}

    # 1. Set the compiler.
    com_opt["CXX"] = os.path.join(paths["MPI_PATH"], "bin", compilers["CXX_MPI"]) if kwargs["mpi"] else compilers["CXX"]

    # 2. Set the OpenMP flags.
    if not kwargs["openmp"]: flags["OPENMPFLAG"] = ""

    # 3. Set the nvcc common flags
    # NOTE: `-G` may cause the GPU Poisson solver to fail
    if kwargs["debug"]: flags["NVCCFLAG_COM"] += "-g -Xptxas -v"

    # 4. Write flags to compile option dictionary.
    for key, val in flags.items():
        com_opt[key] = val

    return com_opt

def validation( paths, depends, constraints, **kwargs ):
    success = True

    # 0. Checking the Makefile_base existence.
    if not os.path.isfile( GAMER_MAKE_BASE ):
        LOGGER.error("%s does not exist."%(GAMER_MAKE_BASE))
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
                LOGGER.error("The option <--%s=%s> requires <--%s> to be set to [%s]. Current: <--%s=%s>."%(opt, str(kwargs[opt]), check_opt, val_str, check_opt, kwargs[check_opt]))
                success = False

    # 2. Checking other conditions.
    # A. Physics
    # A.1 Module
    if kwargs["model"] == "HYDRO":
        if kwargs["passive"] < 0:
            LOGGER.error("Passive scalar should not be negative. Current: %d"%kwargs["passive"])
            success = False

        if kwargs["dual"] not in [NONE_STR, "DE_ENPY"]:
            LOGGER.error("This dual energy form is not supported yet. Current: %s"%kwargs["dual"])
            success = False

    elif kwargs["model"] == "ELBDM":
        if kwargs["passive"] < 0:
            LOGGER.error("Passive scalar should not be negative. Current: %d"%kwargs["passive"])
            success = False

    elif kwargs["model"] == "PAR_ONLY":
        LOGGER.error("<--model=PAR_ONLY> is not supported yet.")
        success = False
    else:
        LOGGER.error("Unrecognized model: %s. Please add to the model choices."%kwargs["model"])
        success = False

    # A.2 Particle
    if kwargs["particle"]:
        if kwargs["star_formation"] and kwargs["store_par_acc"] and not kwargs["store_pot_ghost"]:
            LOGGER.error("<--store_pot_ghost> must be enabled when <--star_formation> and <--store_par_acc> are enabled.")
            success = False
        if not kwargs["gravity"] and not kwargs["tracer"]:
            LOGGER.error("At least one of <--gravity> or <--tracer> must be enabled for <--particle>.")
            success = False
        if kwargs["par_attribute"] < 0:
            LOGGER.error("Number of particle attributes should not be negative. Current: %d"%kwargs["par_attribute"])
            success = False

    # B. Miscellaneous options
    if kwargs["nlevel"] < 1:
        LOGGER.error("<--nlevel> should be greater than zero. Current: %d"%kwargs["nlevel"])
        success = False

    if kwargs["max_patch"] < 1:
        LOGGER.error("<--max_patch> should be greater than zero. Current: %d"%kwargs["max_patch"])
        success = False

    if kwargs["patch_size"]%2 != 0 or kwargs["patch_size"] < 8:
        LOGGER.error("<--patch_size> should be an even number greater than or equal to 8. Current: %d"%kwargs["patch_size"])
        success = False

    if kwargs["overlap_mpi"]:
        LOGGER.error("<--overlap_mpi> is not supported yet.")
        success = False

    if not success: raise BaseException( "The above vaildation failed." )
    return

def warning( paths, **kwargs ):
    # 1. Makefile
    if os.path.isfile( GAMER_MAKE_OUT ):
        LOGGER.warning("%s already exists and will be overwritten."%(GAMER_MAKE_OUT))

    # 2. Physics
    if kwargs["model"] == "ELBDM" and kwargs["passive"] != 0:
        LOGGER.warning("Not supported yet and can only be used as auxiliary fields.")

    # 3. Path
    path_links = { "gpu":{True:"CUDA_PATH"}, "fftw":{"FFTW2":"FFTW2_PATH", "FFTW3":"FFTW3_PATH"},
                   "mpi":{True:"MPI_PATH"}, "hdf5":{True:"HDF5_PATH"}, "grackle":{True:"GRACKLE_PATH"},
                   "gsl":{True:"GSL_PATH"}, "libyt":{True:"LIBYT_PATH"} }

    for arg, links in path_links.items():
        for val, p_name in links.items():
            if kwargs[arg] != val: continue
            if paths.setdefault(p_name, "") != "": continue
            LOGGER.warning("%-15s is not given in %s.config when setting <--%s=%s>"%(p_name, kwargs["machine"], arg, str(val)))
    return



####################################################################################################
# Main execution
####################################################################################################
if __name__ == "__main__":
    # 0. Set the logger
    ch = logging.StreamHandler()
    ch.setFormatter( CustomFormatter() )
    LOGGER.addHandler( ch )

    # 1. Get the execution command
    command = " ".join(["# This makefile is generated by the following command:", "\n#", sys.executable] + sys.argv + ["\n"])
    LOGGER.info( " ".join( [sys.executable] + sys.argv ) )

    # 2. Load the input arguments
    args, name_table, depends, constraints = load_arguments()

    # 3. Prepare the makefile args
    # 3.1 Load the machine setup
    paths, compilers, flags, gpus = load_config( os.path.join(GAMER_CONFIG_DIR, args["machine"]+".config") )

    # 3.2 Validate arguments
    validation( paths, depends, constraints, **args )

    warning( paths, **args )

    # 3.3 Add the SIMU_OPTION
    LOGGER.info("========================================")
    LOGGER.info("GAMER has the following setting.")
    LOGGER.info("----------------------------------------")
    sims = set_sims( name_table, depends, **args )

    # 3.4 Set the compiler
    compiles = set_compile( paths, compilers, flags, args )

    # 3.5 Set the GPU
    gpu_setup = set_gpu( gpus, flags, args )

    # 4. Create Makefile
    # 4.1 Read
    with open( GAMER_MAKE_BASE, "r" ) as make_base:
        makefile = make_base.read()

    # 4.2 Replace
    LOGGER.info("----------------------------------------")
    for key, val in paths.items():
        LOGGER.info("%-25s : %s"%(key, val))
        makefile, num = re.subn(r"@@@%s@@@"%(key), val, makefile)
        if num == 0: raise BaseException("The string @@@%s@@@ is not replaced correctly."%key)

    for key, val in sims.items():
        makefile, num = re.subn(r"@@@%s@@@"%(key), val, makefile)
        if num == 0: raise BaseException("The string @@@%s@@@ is not replaced correctly."%key)

    LOGGER.info("----------------------------------------")
    for key, val in compiles.items():
        LOGGER.info("%-25s : %s"%(key, val))
        makefile, num = re.subn(r"@@@%s@@@"%(key), val, makefile)
        if num == 0: raise BaseException("The string @@@%s@@@ is not replaced correctly."%key)

    LOGGER.info("----------------------------------------")
    for key, val in gpu_setup.items():
        LOGGER.info("%-25s : %s"%(key, val))
        makefile, num = re.subn(r"@@@%s@@@"%(key), val, makefile)
        if num == 0: raise BaseException("The string @@@%s@@@ is not replaced correctly."%key)

    LOGGER.info("----------------------------------------")
    for key in re.findall(r"@@@(.+?)@@@", makefile):
        makefile, num = re.subn(r"@@@%s@@@"%key, "", makefile)
        if num == 0: raise BaseException("The string @@@%s@@@ is not replaced correctly."%key)
        LOGGER.warning("@@@%s@@@ is replaced to '' since there is no given value."%key)

    # 4.3 Write
    with open( GAMER_MAKE_OUT, "w") as make_out:
        make_out.write( command + makefile )

    LOGGER.info("========================================")
    LOGGER.info("%s is created."%GAMER_MAKE_OUT)
    LOGGER.info("========================================")

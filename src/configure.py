"""
Cluster and flags setup
    --cluster            Select the cluster config.  
    --flags              Compiler flags config.

A. Physical models and options of diffierent physical models 
    --model              Select the physical model.
    --passive            Set the number of passive scalars.

A.1 Hydro options
    --flu_scheme         Select the fluid solver for HYDRO model.
    --slope              Select the spatial data reconstruction method.
    --flux               Select the Riemann solver.
    --dual               Select the dual-energy formalism.
    --mhd                Enable magnetohydrodynamic.
    --cosmic_ray         Enable cosmic rays.
    --eos                Select the equation of state.
    --barotropic         Whether or not the --eos set is barotropic.

A.2 ELBDM scheme
    --conserve_mass      Enforce the mass conservation. 
    --laplacian_four     Enable the fourth order of Laplacian.
    --self_interaction   Including the quartic self-interaction potential.

A.3 gravity
    --gravity            Enable gravity.
    --pot_scheme         Select the Poisson solver.
    --store_pot_ghost    Store the potential ghost-zone for each patch on each side.
    --unsplit_gravity    Use unsplitting method to couple gravity to the target model.
    --comoving           Comoving frame for cosmological simulation.

A.4 particle
    --particle           Enable particle.
    --tracer             Enable tracer particles. 
    --store_acc          Store particle acceleration.
    --star_formation     Allow creating new particles after initialization.
    --par_attribute      Set user defined particle attributes.

A.5 grackle
    --grackle            Enable Grackle, a chemistry and radiative cooling library. 

B. Miscellaneous options
    --nlevel             Set the maximum level of AMR.
    --max_patch          Set the maximum patchs on each level of AMR.
    --patch_size         Set size of each direction of a single patch.
    --debug              Enable debug mode.
    --bitwise_reproduce  Enable bitwise reproducibility.
    --timing             Enable to measure timing.
    --timing_solver      Enable measure GPU time.
    --double             Enable double precision.
    --laohu              Work on the NAOC Laohu GPU cluster. 
    --hdf5               Support HDF5 format.
    --GSL                Support GNU scientific library.
    --FFTW               Support FFTW library.
    --LIBYT              Support yt inline analysis.
    --LIBYT_patch        Use patch group as the unit in libyt. Note that this will speed up inline-analysis but increase memory consumption.
    --RNG                Select the random number generator.

C. Parallelization and flags
    --serial             Serial compiler type.
    --openmp             Enable openmp parallization. 
    --mpi                Enable mpi parallization.
    --overlap_mpi        Overlap MPI communication with computation.
    --GPU                Enable GPU.
    --GPU_arch           Select the archtecture of GPU.
"""

####################################################################################################
# Packages
####################################################################################################
import argparse
import re



####################################################################################################
# Global variables
####################################################################################################
GAMER_CONFIG_DIR  = "../configs/"
GAMER_MAKE_BASE   = "Makefile_base"
GAMER_MAKE_OUT    = "Makefile"
GAMER_DESCRIPTION = "Prepare custom Makefile for GAMER."
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
                     "serial":"SERIAL", "openmp":"OPENMP", "mpi":"LOAD_BALANCE=HILBERT", 
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
# Global variables
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

    # C. parallel options
    if kwargs["openmp"]: sim_opt[NAME_TABLE["openmp"]] = kwargs["openmp"]
    if kwargs["mpi"]:
        sim_opt[NAME_TABLE["mpi"]] = kwargs["mpi"]
    else:
        sim_opt["SERIAL"] = True

    if kwargs["GPU"]:
        sim_opt[NAME_TABLE["GPU"]] = kwargs["GPU"]
        sim_opt[NAME_TABLE["GPU_arch"]] = kwargs["GPU_arch"]

    # D. Setup the sumulation option string.
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
            exit("Unknown type to add the simulation options.")

    return {"SIMU_OPTION":opts}

def load_compile( paths, flags, kwargs ):
    com_opt = {}

    # 1. complier. If mpi, overwrite the compiler to mpi.
    com_opt["CXX"] = paths["MPI_PATH"] + "/bin/mpicxx" if kwargs["mpi"] else kwargs["serial"]

    # 2. compiler flags
    if kwargs["serial"] == "icpc" and not kwargs["openmp"]:
        flags["CXXFLAG"] += "-Wno-unknown-pragmas -diag-disable 3180 "

    elif kwargs["serial"] == "g++" and not kwargs["openmp"]:
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
    sim_opt = {}
    
    # A. Physics
    # A.1 Module 
    if   kwargs["model"] == "HYDRO":
        if kwargs["mhd"]: 
            if kwargs["flu_scheme"] not in ["MHM_RP", "CTU"]:
                color_print("MHD only support MHM_RP and CTU.", BCOLOR.FAIL)
                success = False
            if kwargs["flux"] not in ["EXACT", "ROE", "HLLE", "HLLD"]:
                color_print("MHD only support EXACT, ROE, HLLE, and HLLD Riemann solver.", BCOLOR.FAIL)
                success = False
        else:
            if kwargs["flux"] not in ["EXACT", "ROE", "HLLE", "HLLC"]:
                color_print("Pure hydro only support EXACT, ROE, HLLE, and HLLC Riemann solver.", BCOLOR.FAIL)
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
                color_print("ENPY dual energy only support for --eos=GAMMA.", BCOLOR.FAIL)
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
                color_print("ENPY dual energy only support for --cosmic_ray.", BCOLOR.FAIL)
                success = False
        
        if kwargs["barotropic"]:
            if kwargs["eos"] not in ["ISOTHERMAL", "TABULAR", "USER"]:
                color_print("--barotropic is supported for ISOTHERMAL, TABULAR, and USER.", BCOLOR.FAIL)
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
            color_print("Not supported yet and thus can only be used as auxiliary fields.", BCOLOR.FAIL)
            success = False

    elif kwargs["model"] == "PAR_ONLY":
        color_print("--model PAR_ONLY is not supported yet.", BCOLOR.FAIL)
        success = False
    else:
        color_print("Unrecognize model: %s. Please add to the model choices.", BCOLOR.FAIL)
        success = False

    # A.2 Gravity
    if kwargs["gravity"]:
        if not kwargs["FFTW"]:
            color_print("--FFTW must be enable with --gravity.", BCOLOR.FAIL)
            success = False
        if kwargs["unsplit_gravity"] and kwargs["model"] != "HYDRO":
            color_print("--unsplit_gravity is only supported for --model=HYDRO.", BCOLOR.FAIL)
            success = False
    
    #TODO: unsplit gravity and gravity
     
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
    
    if not success: exit(BCOLOR.FAIL+"The above validations are fail."+BCOLOR.ENDC)
    return


def warning( paths, **kwargs ):
    # 1. serial config not match
    if kwargs["serial"] == "icpc" and kwargs["flags"] == "gnu":
        color_print("Warning: The compiler does not match to the default flag config.", BCOLOR.WARNING)
    
    if kwargs["serial"] == "g++" and kwargs["flags"] == "intel":
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
parser = argparse.ArgumentParser( description = GAMER_DESCRIPTION, 
                                  formatter_class = argparse.RawTextHelpFormatter,
                                  epilog = GAMER_EPILOG )

# cluster and flags setup
parser.add_argument( "--cluster", type=str,
                     default="eureka",
                     help="Select the cluster [eureka, YOUR_CLUSTER_NAME]. \nDefault: %(default)s"
                   )

parser.add_argument( "--flags", type=str,
                     default="intel",
                     help="Compiler flags [intel, gnu, YOUR_FLAG_NAME]. \nDefault: %(default)s"
                   )

# A. physical models and options of diffierent physical models 
parser.add_argument( "--model", type=str,
                     default="HYDRO",
                     choices=["HYDRO", "ELBDM", "PAR_ONLY"],
                     help="Select the physical model. \nDefault: %(default)s"
                   )

parser.add_argument( "--passive", type=int,
                     default=0,
                     help="Set the number of passive scalars. \nDefault: %(default)s"
                   )

# A.1 Hydro options
parser.add_argument( "--flu_scheme", type=str,
                     default="CTU",
                     choices=["RTVD", "MHM", "MHM_RP", "CTU"],
                     help="Select the fluid solver for HYDRO model. \nDefault: %(default)s"
                   )

parser.add_argument( "--slope", type=str,
                     default="PPM",
                     choices=["PLM", "PPM"],
                     help="Select the spatial data reconstruction method. \nDefault: %(default)s"
                   )

parser.add_argument( "--flux", type=str,
                     default="ROE",
                     choices=["EXACT", "ROE", "HLLE", "HLLC", "HLLD"],
                     help="Select the Riemann solver. \nDefault: %(default)s"
                   )

parser.add_argument( "--dual", type=str,
                     default="",
                     choices=["", "ENPY", "EINT"],
                     help="Select the dual-energy formalism. \nDefault: %(default)s"
                   )

parser.add_argument( "--mhd",
                     action="store_true",
                     help="Enable magnetohydrodynamic. \nDefault: %(default)s"
                   )

parser.add_argument( "--cosmic_ray",
                     action="store_true",
                     help="Enable cosmic rays. \nDefault: %(default)s"
                   )

parser.add_argument( "--eos", type=str,
                     default="GAMMA",
                     choices=["GAMMA", "ISOTHERMAL", "NUCLEAR", "TABULAR", "USER"],
                     help="Select the equation of state. \nDefault: %(default)s"
                   )

parser.add_argument( "--barotropic",
                     action="store_true",
                     help="Whether or not the --eos set is barotropic. \nDefault: %(default)s"
                   )

# A.2 ELBDM scheme
parser.add_argument( "--conserve_mass",
                     action="store_true",
                     help="Enforce the mass conservation. \nDefault: %(default)s"
                   )

parser.add_argument( "--laplacian_four",
                     action="store_true",
                     help="Enable the fourth order of Laplacian. \nDefault: %(default)s"
                   )

parser.add_argument( "--self_interaction",
                     action="store_true",
                     help="Including the quartic self-interaction potential. \nDefault: %(default)s"
                   )

# A.3 gravity
parser.add_argument( "--gravity",
                     action="store_true",
                     help="Enable gravity. \nDefault: %(default)s"
                   )

parser.add_argument( "--pot_scheme", type=str,
                     default="SOR",
                     choices=["SOR", "MG"],
                     help="Select the Poisson solver. \nDefault: %(default)s"
                   )

parser.add_argument( "--store_pot_ghost",
                     action="store_true",
                     help="Store the potential ghost-zone for each patch on each side. \nDefault: %(default)s"
                   )

parser.add_argument( "--unsplit_gravity",
                     action="store_true",
                     help="Use unsplitting method to couple gravity to the target model. \nDefault: %(default)s"
                   )

parser.add_argument( "--comoving",
                     action="store_true",
                     help="Comoving frame for cosmological simulation. \nDefault: %(default)s"
                   )

# A.4 particle
parser.add_argument( "--particle",
                     action="store_true",
                     help="Enable particle. \nDefault: %(default)s"
                   )
parser.add_argument( "--tracer",
                     action="store_true",
                     help="Enable tracer particles. \nDefault: %(default)s"
                   )

parser.add_argument( "--store_acc",
                     action="store_true",
                     help="Store particle acceleration. \nDefault: %(default)s"
                   )

parser.add_argument( "--star_formation",
                     action="store_true",
                     help="Allow creating new particles after initialization. \nDefault: %(default)s"
                   )

parser.add_argument( "--par_attribute", type=int,
                     default=0,
                     help="Set the number of user defined particle attributes. \nDefault: %(default)s"
                   )

# A.5 grackle
parser.add_argument( "--grackle",
                     action="store_true",
                     help="Enable Grackle, a chemistry and radiative cooling library. \nDefault: %(default)s"
                   )

# B. miscellaneous options
parser.add_argument( "--nlevel", type=int,
                     default=10,
                     help="Set the maximum level of AMR. \nDefault: %(default)s"
                   )

parser.add_argument( "--max_patch", type=int,
                     default=100000,
                     help="Set the maximum patchs on each level of AMR. \nDefault: %(default)s"
                   )

parser.add_argument( "--patch_size", type=int,
                     default=8,
                     help="Set size of each direction of a single patch. \nDefault: %(default)s"
                   )


parser.add_argument( "--debug",
                     action="store_true",
                     help="Enable debug mode. \nDefault: %(default)s"
                   )

parser.add_argument( "--bitwise_reproduce",
                     action="store_true",
                     help="Enable bitwise reproduce. \nDefault: %(default)s"
                   )

parser.add_argument( "--timing",
                     action="store_true",
                     help="Enable to measure timing. \nDefault: %(default)s"
                   )

parser.add_argument( "--timing_solver",
                     action="store_true",
                     help="Enable measure GPU time. \nDefault: %(default)s"
                   )

parser.add_argument( "--double",
                     action="store_true",
                     help="Enable double precision. \nDefault: %(default)s"
                   )

parser.add_argument( "--laohu",
                     action="store_true",
                     help="Work on the NAOC Laohu GPU cluster. \nDefault: %(default)s"
                   )

parser.add_argument( "--hdf5",
                     action="store_true",
                     help="Support HDF5 format. \nDefault: %(default)s"
                   )

parser.add_argument( "--GSL",
                     action="store_true",
                     help="Support GNU scientific library. \nDefault: %(default)s"
                   )

parser.add_argument( "--FFTW",
                     action="store_true",
                     help="Support FFTW library. \nDefault: %(default)s"
                   )

parser.add_argument( "--LIBYT",
                     action="store_true",
                     help="Support yt inline analysis. \nDefault: %(default)s"
                   )

parser.add_argument( "--LIBYT_patch",
                     action="store_true",
                     help="Use patch group as the unit in libyt. Note that this will speed up inline-analysis but increase memory consumption. \nDefault: %(default)s"
                   )

parser.add_argument( "--RNG", type=str,
                     default="RNG_GNU_EXT",
                     choices=["RNG_GNU_EXT", "RNG_CPP11"],
                     help="Select the random number generator. \nDefault: %(default)s"
                   )

# C. parallelization and flags
parser.add_argument( "--serial",
                     default="icpc",
                     help="Serial compiler type[icpc, g++, YOUR_GNU_PATH]. \nDefault: %(default)s"
                   )
parser.add_argument( "--openmp",
                     action="store_true",
                     help="Enable openmp parallization. \nDefault: %(default)s"
                   )
 
parser.add_argument( "--mpi",
                     action="store_true",
                     help="Enable mpi parallization. \nDefault: %(default)s"
                   )

parser.add_argument( "--overlap_mpi",
                     action="store_true",
                     help="Overlap MPI communication with computation. \nDefault: %(default)s"
                   )

parser.add_argument( "--GPU",
                     action="store_true",
                     help="Enable GPU. \nDefault: %(default)s"
                   )

parser.add_argument( "--GPU_arch", type=str,
                     default="TURING",
                     choices=["FERMI", "KEPLER", "MAXWELL", "PASCAL", "VOLTA", "TURING", "AMPERE"],
                     help="Select the archtecture of GPU. \nDefault: %(default)s"
                   )


args = vars( parser.parse_args() )


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

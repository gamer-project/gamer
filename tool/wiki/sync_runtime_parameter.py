#!/bin/python3
import re
from string import ascii_uppercase as auc



#====================================================================================================
# Global constants
#====================================================================================================
ALL_PARAM_FILE = "../../example/input/Input__Parameter"
PARAM_CPP_FILE = "../../src/Init/Init_Load_Parameter.cpp"
LINK_FILES     = [ "../../doc/wiki/Runtime-Parameters:-General.md",
                   "../../doc/wiki/MPI-and-OpenMP.md",
                   "../../doc/wiki/GPU.md",
                   "../../doc/wiki/Runtime-Parameters:-Units.md",
                   "../../doc/wiki/Initial-Conditions.md",
                   "../../doc/wiki/Hydro.md",
                   "../../doc/wiki/Gravity.md",
                   "../../doc/wiki/Particles.md",
                   "../../doc/wiki/Runtime-Parameters:-Cosmology.md",
                   "../../doc/wiki/Chemistry-and-Radiation.md",
                   "../../doc/wiki/Star-Formation.md",
                   "../../doc/wiki/Feedback.md",
                   "../../doc/wiki/Runtime-Parameters:-Timestep.md",
                   "../../doc/wiki/Runtime-Parameters:-Refinement.md",
                   "../../doc/wiki/Runtime-Parameters:-Interpolation.md",
                   "../../doc/wiki/Outputs.md",
                   "../../doc/wiki/Runtime-Parameters:-Miscellaneous.md"
                 ]
OUT_MD         = "Runtime-Parameters:-All.md"
REPLACE_DICT   = { "NoMin_double":"None", "NoMax_double":"None", "NoDef_double":"None",
                   "NoMin_int":"None", "NoMax_int":"None",
                   "NoMin_long":"None", "NoMax_long":"None",
                   "NoDef_str":"None", "Useless_str":"None",
                   "Useless_bool":"None",
                   "true":"1", "false":"0",
                   "__DBL_MAX__":"1.79769313e+308", "Eps_double":"2.22507386e-308"
                 }



#====================================================================================================
# Classes
#====================================================================================================
class parameter():
    def __init__( self, string ):
        name, description = self.get_name_description( string )
        self.name         = name
        self.link_name    = name
        self.default      = []
        self.minimum      = []
        self.maximum      = []
        self.description  = description
        self.NAdd         = 0

    def get_name_description( self, string ):
        words       = string.split()
        name        = words[0]
        description = ""
        for i in range(3, len(words)):
            description += words[i] + ' '
        return name, description[:-1]

    def append_description( self, string ):
        words = string.split()
        self.description += "<br />" if self.name == "TESTPROB_ID" else ' '
        self.description += ' '.join(words[1:])
        return

    def get_link_name( self, file_dict ):
        for file_name in file_dict:
            search_pattern = "[%s](#%s)"%(self.name, self.name)
            if search_pattern in file_dict[file_name]:
                self.link_name = "[[ %s \\| %s#%s ]]"%(self.name, file_name, self.name)
                return True
        return False



#====================================================================================================
# Main
#====================================================================================================
with open( ALL_PARAM_FILE, 'r' ) as f:
    lines = f.readlines()

with open( PARAM_CPP_FILE, 'r' ) as f:
    lines_cpp = f.readlines()

link_source_mds = {}
for file_path in LINK_FILES:
    file_name = file_path.split('/')[-1][:-3] # get the file name without path and trailing `.md`
    with open( file_path, 'r' ) as f:
        link_source_mds[file_name] = f.read()

params = {}

# get all parameters from ALL_PARAM_FILE
# NOTE: assuming the line starts with space should be the comment of the above line
last_key = ""
for i, line in enumerate(lines):
    if line[0] == '\n': continue
    if line[0] == '#': continue
    if line[0] == ' ':
        params[last_key].append_description( line )
        continue
    last_key = line.split()[0]
    params[last_key] = parameter( line )

# get all default, min, max, value from PARAM_CPP_FILE
for i, line in enumerate(lines_cpp):
    if "ReadPara->Add" not in line: continue
    words = list( filter( None, re.split( ',| ', line ) ) )
    if "//" == words[0]:  continue
    key     = words[1][1:-1]
    default = words[3]
    minimum = words[4]
    maximum = words[5]
    try:
        if params[key].NAdd >= 1: print( "%s has two or more Add function"%key )
        params[key].default.append( REPLACE_DICT[default] if default in REPLACE_DICT else default )
        params[key].minimum.append( REPLACE_DICT[minimum] if minimum in REPLACE_DICT else minimum )
        params[key].maximum.append( REPLACE_DICT[maximum] if maximum in REPLACE_DICT else maximum )
        params[key].NAdd += 1
    except:
        print( key, "does not exist in %s"%ALL_PARAM_FILE )

# get the detailed description link from LINK_FILES
for p in params:
    status = params[p].get_link_name( link_source_mds )

# output markdown file
with open( OUT_MD, 'w' ) as f:
    param_str_format = '| %-100s | %15s | %15s | %15s | %s |\n'

    f.write( '> [!CAUTION]\n' )
    f.write( '> Please do not edit this file(page) manually since the workflow will overwrite your changes.\n' )
    f.write( '\n' )
    f.write( 'This file(page) is automatically generated by the workflow `Update all parameters wiki page` using the script `tool/wiki/sync_runtime_parameter.py`.\n' )
    f.write( '\n' )
    f.write( 'The workflow is triggered by push changes to any of  `src/Init/Init_Load_Parameter.cpp`, `example/input/Input__Paramter`, and `.github/workflow/sync_runtime_parameter.py`.\n' )
    f.write( '\n' )
    f.write( '# Index\n' )
    f.write( ', '.join( ['[%s](#%s)'%(i, i) for i in auc] ) )
    f.write( '\n' )

    start_char = ''
    params_sorted_key = sorted( params.keys() )
    for key in params_sorted_key:
        # Add alphabet title
        if start_char != key[0]:
            start_char = key[0]
            f.write( '\n' )
            f.write( "# %s\n"%(start_char) )
            f.write( param_str_format%("Name", "Default", "Min", "Max", "Short description") )
            f.write( param_str_format%(":---", ":---", ":---", ":---", ":---") )

        for i in range(params[key].NAdd):
            string = param_str_format%(params[key].link_name, params[key].default[i], params[key].minimum[i], params[key].maximum[i], params[key].description)
            f.write( string )

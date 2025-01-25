#!/bin/python3
"""
## What it can do

Automatically update `doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-All.md`
when any of the following files are changed, updated, or pushed:
- `tool/wiki/sync_runtime_parameter.py`
- `src/Init/Init_Load_Parameter.cpp`
- `example/test_problem/Template/Input__Parameter`

## How it works

1. Detection:
   The new workflow, "Update All Parameters Wiki Page", detects changes in any of the three specified files.
2. Execution:
   If changes are detected, the workflow runs the script `tool/wiki/sync_runtime_parameter.py`.
3. Update:
   The script updates `doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-All.md` by creating
   a new commit with the message: `[Workflow] Update all parameters wiki page`.

## Script Algorithm

1. Extract Parameter Names and Comments:
   Retrieve all parameter names and their associated short comments from `example/test_problem/Template/Input__Parameter`.
2. Retrieve Parameter Details:
   Gather the minimum, maximum, and default values of each parameter from `src/Init/Init_Load_Parameter.cpp`.
3. Link Integration:
   Identify the original page links using the references placed at the top of each runtime parameter
   wiki page (e.g., `[BOX_SIZE](#BOX_SIZE), &nbsp;`).
4. Markdown Generation:
   Compile the extracted data and generate the updated markdown file.
"""
#====================================================================================================
# Imports
#====================================================================================================
import re
from string import ascii_uppercase as auc



#====================================================================================================
# Global constants
#====================================================================================================
ALL_PARAM_FILE = "../../example/input/Input__Parameter"
PARAM_CPP_FILE = "../../src/Init/Init_Load_Parameter.cpp"
LINK_FILES     = [ "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Chemistry-and-Radiation.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Cosmology.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Feedback.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-GPU.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-General.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Gravity.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Hydro.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Initial-Conditions.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Interpolation.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-MPI-and-OpenMP.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Miscellaneous.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Outputs.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Particles.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Refinement.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Star-Formation.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Timestep.md",
                   "../../doc/wiki/Runtime-Parameters-related/Runtime-Parameters:-Units.md"
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
        self.default      = ""
        self.minimum      = ""
        self.maximum      = ""
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
    file_name = file_path.split('/')[-1][:-3] # get the filename without path and trailing `.md`
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

# get all default, min, max, and value from PARAM_CPP_FILE
for i, line in enumerate(lines_cpp):
    if "ReadPara->Add" not in line: continue
    words = list( filter( None, re.split( ',| ', line ) ) )
    if "//" == words[0]:  continue
    key     = words[1][1:-1]
    default = words[3]
    minimum = words[4]
    maximum = words[5]

    if key not in params:
        print( "%-30s does not exist in %s"%(key, ALL_PARAM_FILE) )
        continue

    if params[key].NAdd >= 1:
        print( "%-30s has more than one ReadPara->Add() function in %s"%(key, PARAM_CPP_FILE) )
        params[key].default = "Depend"
        params[key].minimum = "Depend"
        params[key].maximum = "Depend"
    else:
        params[key].default = REPLACE_DICT[default] if default in REPLACE_DICT else default
        params[key].minimum = REPLACE_DICT[minimum] if minimum in REPLACE_DICT else minimum
        params[key].maximum = REPLACE_DICT[maximum] if maximum in REPLACE_DICT else maximum
    params[key].NAdd += 1

# get the detailed description link from LINK_FILES
for p in params:
    status = params[p].get_link_name( link_source_mds )
    if not status:
        print( "Can not find the description of %-30s in all runtime parameter pages!"%(p) )

# output the markdown file
with open( OUT_MD, 'w' ) as f:
    param_str_format = '| %-100s | %15s | %15s | %15s | %s |\n'

    f.write( '> [!CAUTION]\n' )
    f.write( '> Please do not edit this file (page) manually, since the workflow will overwrite your changes.\n' )
    f.write( '\n' )
    f.write( 'This file (page) is automatically generated by the workflow `Update all parameters wiki page` using the script `tool/wiki/sync_runtime_parameter.py`.\n' )
    f.write( '\n' )
    f.write( 'The workflow is triggered when changes are pushed to any of the following files:\n' )
    f.write( '- `src/Init/Init_Load_Parameter.cpp`\n' )
    f.write( '- `example/input/Input__Paramter`\n' )
    f.write( '- `tool/wiki/sync_runtime_parameter.py`\n' )
    f.write( '\n' )
    f.write( 'For variables with `Default/Min/Max` labeled as `Depend`, click the parameter names for more details.\n' )
    f.write( '\n' )
    f.write( '# Index\n' )
    f.write( ', '.join( ['[%s](#%s)'%(i, i) for i in auc] ) )
    f.write( '\n' )

    start_char = ''
    params_sorted_key = sorted( params.keys() )
    for key in params_sorted_key:
        # add alphabet title
        if start_char != key[0]:
            start_char = key[0]
            f.write( '\n' )
            f.write( "# %s\n"%(start_char) )
            f.write( param_str_format%("Name", "Default", "Min", "Max", "Short description") )
            f.write( param_str_format%(":---", ":---", ":---", ":---", ":---") )

        string = param_str_format%(params[key].link_name, params[key].default, params[key].minimum, params[key].maximum, params[key].description)
        f.write( string )

    f.write( '\n' )
    f.write( '\n' )
    f.write( '## Remarks\n' )
    f.write( '\n' )
    f.write( '\n' )
    f.write( '<br>\n' )
    f.write( '\n' )
    f.write( '## Links\n' )
    f.write( '* [[Main page of Runtime Parameters | Runtime Parameters]]\n' )

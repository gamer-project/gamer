import argparse
import sys
import numpy as np


# load the command-line parameters
parser = argparse.ArgumentParser( description='Analyze the average evolution time-step at the specified AMR level' )

parser.add_argument( '-l', action='store', required=True,  type=int, dest='lv',
                     help='AMR level' )
parser.add_argument( '-n', action='store', required=False, type=int, dest='nave',
                     help='number of sub-steps to average over [%(default)d]', default=10 )
parser.add_argument( '-c', action='store', required=False, type=int, dest='column',
                     help='target time-step column [%(default)d]', default=5 )
parser.add_argument( '-i', action='store', required=False, type=str, dest='filename_in',
                     help='filename of the simulation time-step log file [%(default)s]', default='Record__TimeStep' )
parser.add_argument( '-o', action='store', required=True,  type=str, dest='filename_out',
                     help='output filename' )

args=parser.parse_args()

# check
assert args.lv     >= 0, '-l (%d) < 0' % (args.lv)
assert args.column >= 0, '-c (%d) < 0' % (args.column)
assert args.nave   >= 1, '-n (%d) < 1' % (args.nave)


# take note
File_Out = open( args.filename_out, "w" )

File_Out.write( '#Command-line arguments:\n' )
File_Out.write( '#-------------------------------------------------------------------\n' )
File_Out.write( '#' )
for t in range( len(sys.argv) ):
   File_Out.write( ' %s' % str(sys.argv[t]) )
File_Out.write( '\n' )
File_Out.write( '#-------------------------------------------------------------------\n\n' )



# load the level, time, and dt from the simulation log file
log = np.loadtxt( args.filename_in, usecols=(0,3,args.column), unpack=True )


# get the time and dt at the target level
row_lv = log[ 0, : ] == args.lv
t      = log[ 1, row_lv ]
dt     = log[ 2, row_lv ]


# calculate and record the average time and dt
t_ave  = 0.0
dt_ave = 0.0

File_Out.write( "#%13s   %13s\n" % ("Time", "dt") )

for i in range( 0, dt.size, 1 ):
   t_ave  += t [i]
   dt_ave += dt[i]

   if (i+1)%args.nave == 0:
      t_ave  = t_ave  / args.nave
      dt_ave = dt_ave / args.nave

      File_Out.write( "% 13.7e   %13.7e\n" % (t_ave, dt_ave) )

      t_ave  = 0.0
      dt_ave = 0.0

File_Out.close()


import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Output the center' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='data path prefix [%(default)s]', default='../' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix


def _TotDens(field, data):
    return data["Dens"]+data["ParDens"]

yt.add_field( ("gamer","TotDens"),
              function=_TotDens, units="code_mass/code_length**3",
              sampling_type="cell")


field_dens       = ('gamer','Dens')
field_pardens    = ('gamer','ParDens')
field_totdens    = ('gamer','TotDens')
field_pote       = ('gamer','Pote')


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
    print( '-------------------------------------------------------------------' )
    print("")
    print("Data name              :", ds )
    print("Time                   :", ds.current_time.in_units("code_time").d )
    print("")

    ad = ds.all_data()

    #Extrema_density     = ad.quantities.extrema(field_dens).in_units("code_density")
    #print("Extrema_density:", Extrema_density )
    MaxLoc_density      = ad.quantities.max_location(field_dens)
    print("Max Density Value      :", MaxLoc_density[0].in_units("code_density").d )
    print("Max Density Coord_x    :", MaxLoc_density[1].in_units("code_length").d )
    print("Max Density Coord_y    :", MaxLoc_density[2].in_units("code_length").d )
    print("Max Density Coord_z    :", MaxLoc_density[3].in_units("code_length").d )
    print("")

    #Extrema_pardensity  = ad.quantities.extrema(field_pardens).in_units("code_density")
    #print("Extrema_pardensity:", Extrema_pardensity )
    MaxLoc_pardensity   = ad.quantities.max_location(field_pardens)
    print("Max Par_Density Value  :", MaxLoc_pardensity[0].in_units("code_density").d )
    print("Max Par_Density Coord_x:", MaxLoc_pardensity[1].in_units("code_length").d )
    print("Max Par_Density Coord_y:", MaxLoc_pardensity[2].in_units("code_length").d )
    print("Max Par_Density Coord_z:", MaxLoc_pardensity[3].in_units("code_length").d )
    print("")

    #Extrema_totdensity  = ad.quantities.extrema(field_totdens).in_units("code_density")
    #print("Extrema_totdensity:", Extrema_totdensity )
    MaxLoc_totdensity   = ad.quantities.max_location(field_totdens)
    print("Max Tot_Density Value  :", MaxLoc_totdensity[0].in_units("code_density").d )
    print("Max Tot_Density Coord_x:", MaxLoc_totdensity[1].in_units("code_length").d )
    print("Max Tot_Density Coord_y:", MaxLoc_totdensity[2].in_units("code_length").d )
    print("Max Tot_Density Coord_z:", MaxLoc_totdensity[3].in_units("code_length").d )
    print("")

    #Extrema_pote        = ad.quantities.extrema(field_pote)
    #print("Extrema_pote:", Extrema_pote )
    #MinLoc_pote         = ad.quantities.min_location(field_pote)
    #print("Min Potential Value    :", MinLoc_pote[0].in_units("code_length**2/code_time**2").d )
    #print("Min Potential Coord_x  :", MinLoc_pote[1].in_units("code_length").d )
    #print("Min Potential Coord_y  :", MinLoc_pote[2].in_units("code_length").d )
    #print("Min Potential Coord_z  :", MinLoc_pote[3].in_units("code_length").d )
    #print("")

    CoM                 = ad.quantities.center_of_mass(use_gas=True, use_particles=True, particle_type='all')
    print("CoM Coord_x            :", CoM[0].in_units("code_length").d )
    print("CoM Coord_y            :", CoM[1].in_units("code_length").d )
    print("CoM Coord_z            :", CoM[2].in_units("code_length").d )
    print("")
    print( '-------------------------------------------------------------------' )


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
    return data['Dens']+data['ParDens']


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():
    print( '' )
    print( '-------------------------------------------------------------------' )
    print( 'Data name               = ', ds )
    print( 'Time                    = % 14.7e'%( ds.current_time.in_units('code_time').d ) )
    print( '-------------------------------------------------------------------' )
    print( '' )

    if ds.parameters["Particle"] == 1:
        ds.add_field( ('gamer', 'TotDens'), function=_TotDens, units='code_mass/code_length**3', sampling_type='cell' )

    ad = ds.all_data()

    if ds.parameters["Model"] == 'Hydro' or ds.parameters["Model"] == 'ELBDM':
        #Extrema_density     = ad.quantities.extrema( ('gamer','Dens') ).in_units('code_density')
        #print( 'Extrema_density         = ', Extrema_density )
        MaxLoc_density      = ad.quantities.max_location( ('gamer','Dens') )
        print( 'Max Density Value       = % 14.7e'%( MaxLoc_density[0].in_units('code_density').d ) )
        print( 'Max Density Coord_x     = % 14.7e'%( MaxLoc_density[1].in_units('code_length').d  ) )
        print( 'Max Density Coord_y     = % 14.7e'%( MaxLoc_density[2].in_units('code_length').d  ) )
        print( 'Max Density Coord_z     = % 14.7e'%( MaxLoc_density[3].in_units('code_length').d  ) )
        print( '' )

    if ds.parameters["Particle"] == 1:
        #Extrema_pardensity  = ad.quantities.extrema( ('gamer','ParDens') ).in_units('code_density')
        #print( 'Extrema_pardensity      = ', Extrema_pardensity )
        MaxLoc_pardensity   = ad.quantities.max_location( ('gamer','ParDens') )
        print( 'Max Par_Density Value   = % 14.7e'%( MaxLoc_pardensity[0].in_units('code_density').d ) )
        print( 'Max Par_Density Coord_x = % 14.7e'%( MaxLoc_pardensity[1].in_units('code_length').d  ) )
        print( 'Max Par_Density Coord_y = % 14.7e'%( MaxLoc_pardensity[2].in_units('code_length').d  ) )
        print( 'Max Par_Density Coord_z = % 14.7e'%( MaxLoc_pardensity[3].in_units('code_length').d  ) )
        print( '' )

        #Extrema_totdensity  = ad.quantities.extrema( ('gamer','TotDens') ).in_units('code_density')
        #print( 'Extrema_totdensity      = ', Extrema_totdensity )
        MaxLoc_totdensity   = ad.quantities.max_location( ('gamer','TotDens') )
        print( 'Max Tot_Density Value   = % 14.7e'%( MaxLoc_totdensity[0].in_units('code_density').d ) )
        print( 'Max Tot_Density Coord_x = % 14.7e'%( MaxLoc_totdensity[1].in_units('code_length').d  ) )
        print( 'Max Tot_Density Coord_y = % 14.7e'%( MaxLoc_totdensity[2].in_units('code_length').d  ) )
        print( 'Max Tot_Density Coord_z = % 14.7e'%( MaxLoc_totdensity[3].in_units('code_length').d  ) )
        print( '' )

    if ds.parameters["Opt__Output_Pot"] == 1:
        #Extrema_pote        = ad.quantities.extrema( ('gamer','Pote') )
        #print( 'Extrema_pote            = ', Extrema_pote )
        MinLoc_pote         = ad.quantities.min_location( ('gamer','Pote') )
        print( 'Min Potential Value     = % 14.7e'%( MinLoc_pote[0].in_units('code_length**2/code_time**2').d ) )
        print( 'Min Potential Coord_x   = % 14.7e'%( MinLoc_pote[1].in_units('code_length').d                 ) )
        print( 'Min Potential Coord_y   = % 14.7e'%( MinLoc_pote[2].in_units('code_length').d                 ) )
        print( 'Min Potential Coord_z   = % 14.7e'%( MinLoc_pote[3].in_units('code_length').d                 ) )
        print( '' )
    else:
        print( 'WARNING : To find the minimum gravitational potential, please turn on OPT__OUTPUT_POT in Input__Parameter !!\n' )

    if ds.parameters["Model"] == 'Hydro' or ds.parameters["Model"] == 'ELBDM':
        CoM_Gas             = ad.quantities.center_of_mass( use_gas=True, use_particles=False )
        print( 'CoM_Gas Coord_x         = % 14.7e'%( CoM_Gas[0].in_units('code_length').d ) )
        print( 'CoM_Gas Coord_y         = % 14.7e'%( CoM_Gas[1].in_units('code_length').d ) )
        print( 'CoM_Gas Coord_z         = % 14.7e'%( CoM_Gas[2].in_units('code_length').d ) )
        print( '' )

    if ds.parameters["Particle"] == 1:
        CoM_Par             = ad.quantities.center_of_mass( use_gas=False, use_particles=True, particle_type='all' )
        print( 'CoM_Par Coord_x         = % 14.7e'%( CoM_Par[0].in_units('code_length').d ) )
        print( 'CoM_Par Coord_y         = % 14.7e'%( CoM_Par[1].in_units('code_length').d ) )
        print( 'CoM_Par Coord_z         = % 14.7e'%( CoM_Par[2].in_units('code_length').d ) )
        print( '' )

        CoM_All             = ad.quantities.center_of_mass( use_gas=True, use_particles=True, particle_type='all' )
        print( 'CoM_All Coord_x         = % 14.7e'%( CoM_All[0].in_units('code_length').d ) )
        print( 'CoM_All Coord_y         = % 14.7e'%( CoM_All[1].in_units('code_length').d ) )
        print( 'CoM_All Coord_z         = % 14.7e'%( CoM_All[2].in_units('code_length').d ) )
        print( '')


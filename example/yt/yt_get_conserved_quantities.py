import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Output the conserved quantities' )

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


# define new fields
def _kinetic_energy(field, data):
    return data['kinetic_energy_density']*data['cell_volume']

def _potential_energy(field, data):
    return 0.5*data['Pote']*data['Dens']*data['cell_volume'] if data.ds.parameters["Opt__Output_Pot"] == 1 else float('nan')*data.ds.units.code_mass*data.ds.units.code_velocity**2

def _thermal_energy(field, data):
    return data['thermal_energy_density']*data['cell_volume']

def _total_energy(field, data):
    return data['total_energy_density']*data['cell_volume'] + data['potential_energy']

def _particle_momentum_x(field, data):
    return data['particle_mass']*data['particle_velocity_x']

def _particle_momentum_y(field, data):
    return data['particle_mass']*data['particle_velocity_y']

def _particle_momentum_z(field, data):
    return data['particle_mass']*data['particle_velocity_z']

def _particle_kinetic_energy(field, data):
    return 0.5*data['particle_mass']*data['particle_velocity_magnitude']**2

def _particle_potential_energy(field, data):
    return 0.5*data['Pote']*data['ParDens']*data['cell_volume'] if data.ds.parameters["Opt__Output_Pot"] == 1 else float('nan')*data.ds.units.code_mass*data.ds.units.code_velocity**2


# add new fields
yt.add_field( ('gamer','kinetic_energy'),            function=_kinetic_energy,            units='code_mass*code_velocity**2', sampling_type='cell'     )
yt.add_field( ('gamer','potential_energy'),          function=_potential_energy,          units='code_mass*code_velocity**2', sampling_type='cell'     )
yt.add_field( ('gamer','thermal_energy'),            function=_thermal_energy,            units='code_mass*code_velocity**2', sampling_type='cell'     )
yt.add_field( ('gamer','total_energy'),              function=_total_energy,              units='code_mass*code_velocity**2', sampling_type='cell'     )
yt.add_field( ('gamer','particle_momentum_x'),       function=_particle_momentum_x,       units='code_mass*code_velocity',    sampling_type='particle' )
yt.add_field( ('gamer','particle_momentum_y'),       function=_particle_momentum_y,       units='code_mass*code_velocity',    sampling_type='particle' )
yt.add_field( ('gamer','particle_momentum_z'),       function=_particle_momentum_z,       units='code_mass*code_velocity',    sampling_type='particle' )
yt.add_field( ('gamer','particle_kinetic_energy'),   function=_particle_kinetic_energy,   units='code_mass*code_velocity**2', sampling_type='particle' )
yt.add_field( ('gamer','particle_potential_energy'), function=_particle_potential_energy, units='code_mass*code_velocity**2', sampling_type='cell'     )


#yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

    print( '')
    print( '-------------------------------------------------------------------' )
    print( 'Data name   = ', ds )
    print( 'Time        = % 14.7e'%( ds.current_time.in_units('code_time').d ) )
    print( '-------------------------------------------------------------------' )

    ad = ds.all_data()

    # Check if the potential field exist
    if ds.parameters["Opt__Output_Pot"] == 0:
        print( 'WARNING : To calculate the gravitational potential energy, please turn on OPT__OUTPUT_POT in Input__Parameter !!\n' )

    # Gas
    Mass_Gas    = ad.quantities.total_quantity( 'cell_mass'                   ) # total HYDRO mass
    MomX_Gas    = ad.quantities.total_quantity( 'momentum_x'                  ) # total HYDRO momentum x
    MomY_Gas    = ad.quantities.total_quantity( 'momentum_y'                  ) # total HYDRO momentum y
    MomZ_Gas    = ad.quantities.total_quantity( 'momentum_z'                  ) # total HYDRO momentum z
    AngMomX_Gas = ad.quantities.total_quantity( 'angular_momentum_x'          ) # total HYDRO angular momentum x
    AngMomY_Gas = ad.quantities.total_quantity( 'angular_momentum_y'          ) # total HYDRO angular momentum y
    AngMomZ_Gas = ad.quantities.total_quantity( 'angular_momentum_z'          ) # total HYDRO angular momentum z
    Ekin_Gas    = ad.quantities.total_quantity( 'kinetic_energy'              ) # total HYDRO kinetic energy
    Eint_Gas    = ad.quantities.total_quantity( 'thermal_energy'              ) # total HYDRO internal energy
    Epot_Gas    = ad.quantities.total_quantity( 'potential_energy'            ) # total HYDRO potential energy
    Etot_Gas    = ad.quantities.total_quantity( 'total_energy'                ) # total HYDRO energy

    # Particle
    Mass_Par    = ad.quantities.total_quantity( 'particle_mass'               ) # total PARTICLE mass
    MomX_Par    = ad.quantities.total_quantity( 'particle_momentum_x'         ) # total PARTICLE momentum x
    MomY_Par    = ad.quantities.total_quantity( 'particle_momentum_y'         ) # total PARTICLE momentum y
    MomZ_Par    = ad.quantities.total_quantity( 'particle_momentum_z'         ) # total PARTICLE momentum z
    AngMomX_Par = ad.quantities.total_quantity( 'particle_angular_momentum_x' ) # total PARTICLE angular momentum x
    AngMomY_Par = ad.quantities.total_quantity( 'particle_angular_momentum_y' ) # total PARTICLE angular momentum y
    AngMomZ_Par = ad.quantities.total_quantity( 'particle_angular_momentum_z' ) # total PARTICLE angular momentum z
    Ekin_Par    = ad.quantities.total_quantity( 'particle_kinetic_energy'     ) # total PARTICLE kinetic energy
    Epot_Par    = ad.quantities.total_quantity( 'particle_potential_energy'   ) # total PARTICLE potential energy
    Etot_Par    = Ekin_Par + Epot_Par                                           # total PARTICLE energy

    # All
    Mass_All    =    Mass_Gas +    Mass_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE mass
    MomX_All    =    MomX_Gas +    MomX_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE momentum x
    MomY_All    =    MomY_Gas +    MomY_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE momentum y
    MomZ_All    =    MomZ_Gas +    MomZ_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE momentum z
    AngMomX_All = AngMomX_Gas + AngMomX_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE angular momentum x
    AngMomY_All = AngMomY_Gas + AngMomY_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE angular momentum y
    AngMomZ_All = AngMomZ_Gas + AngMomZ_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE angular momentum z
    Ekin_All    =    Ekin_Gas +    Ekin_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE kinetic energy
    Epot_All    =    Epot_Gas +    Epot_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE potential energy
    Etot_All    =    Etot_Gas +    Etot_Par                                     # sum of the total HYDRO/ELBDM + PARTICLE energy


    # print the conserved quantities
    print( '' )
    print( 'Gas' )
    print( 'Mass_Gas    = % 14.7e'%(    Mass_Gas.in_units('code_mass').d                           ) )
    print( 'MomX_Gas    = % 14.7e'%(    MomX_Gas.in_units('code_mass*code_velocity').d             ) )
    print( 'MomY_Gas    = % 14.7e'%(    MomY_Gas.in_units('code_mass*code_velocity').d             ) )
    print( 'MomZ_Gas    = % 14.7e'%(    MomZ_Gas.in_units('code_mass*code_velocity').d             ) )
    print( 'AngMomX_Gas = % 14.7e'%( AngMomX_Gas.in_units('code_length*code_mass*code_velocity').d ) )
    print( 'AngMomY_Gas = % 14.7e'%( AngMomY_Gas.in_units('code_length*code_mass*code_velocity').d ) )
    print( 'AngMomZ_Gas = % 14.7e'%( AngMomZ_Gas.in_units('code_length*code_mass*code_velocity').d ) )
    print( 'Ekin_Gas    = % 14.7e'%(    Ekin_Gas.in_units('code_mass*code_velocity**2').d          ) )
    print( 'Eint_Gas    = % 14.7e'%(    Eint_Gas.in_units('code_mass*code_velocity**2').d          ) )
    print( 'Epot_Gas    = % 14.7e'%(    Epot_Gas.in_units('code_mass*code_velocity**2').d          ) )
    print( 'Etot_Gas    = % 14.7e'%(    Etot_Gas.in_units('code_mass*code_velocity**2').d          ) )

    print( '' )
    print( 'Par' )
    print( 'Mass_Par    = % 14.7e'%(    Mass_Par.in_units('code_mass').d                           ) )
    print( 'MomX_Par    = % 14.7e'%(    MomX_Par.in_units('code_mass*code_velocity').d             ) )
    print( 'MomY_Par    = % 14.7e'%(    MomY_Par.in_units('code_mass*code_velocity').d             ) )
    print( 'MomZ_Par    = % 14.7e'%(    MomZ_Par.in_units('code_mass*code_velocity').d             ) )
    print( 'AngMomX_Par = % 14.7e'%( AngMomX_Par.in_units('code_length*code_mass*code_velocity').d ) )
    print( 'AngMomY_Par = % 14.7e'%( AngMomY_Par.in_units('code_length*code_mass*code_velocity').d ) )
    print( 'AngMomZ_Par = % 14.7e'%( AngMomZ_Par.in_units('code_length*code_mass*code_velocity').d ) )
    print( 'Ekin_Par    = % 14.7e'%(    Ekin_Par.in_units('code_mass*code_velocity**2').d          ) )
    print( 'Epot_Par    = % 14.7e'%(    Epot_Par.in_units('code_mass*code_velocity**2').d          ) )
    print( 'Etot_Par    = % 14.7e'%(    Etot_Par.in_units('code_mass*code_velocity**2').d          ) )

    print( '' )
    print( 'All' )
    print( 'Mass_All    = % 14.7e'%(    Mass_All.in_units('code_mass').d                           ) )
    print( 'MomX_All    = % 14.7e'%(    MomX_All.in_units('code_mass*code_velocity').d             ) )
    print( 'MomY_All    = % 14.7e'%(    MomY_All.in_units('code_mass*code_velocity').d             ) )
    print( 'MomZ_All    = % 14.7e'%(    MomZ_All.in_units('code_mass*code_velocity').d             ) )
    print( 'AngMomX_All = % 14.7e'%( AngMomX_All.in_units('code_length*code_mass*code_velocity').d ) )
    print( 'AngMomY_All = % 14.7e'%( AngMomY_All.in_units('code_length*code_mass*code_velocity').d ) )
    print( 'AngMomZ_All = % 14.7e'%( AngMomZ_All.in_units('code_length*code_mass*code_velocity').d ) )
    print( 'Ekin_All    = % 14.7e'%(    Ekin_All.in_units('code_mass*code_velocity**2').d          ) )
    print( 'Epot_All    = % 14.7e'%(    Epot_All.in_units('code_mass*code_velocity**2').d          ) )
    print( 'Etot_All    = % 14.7e'%(    Etot_All.in_units('code_mass*code_velocity**2').d          ) )

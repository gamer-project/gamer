import yt
import gc



class RevisedDatasets(list):
    def __init__(self, filenames):
        super().__init__(filenames)
        self.filenames = filenames

    def __getitem__(self, i):

        if isinstance(i, slice):
            return RevisedDatasets(self.filenames[i])

        ds0 = yt.load(self.filenames[i])

        # override the comoving units in GAMER, so Mpc/h is physical Mpc/h while Mpccm/h is comoving Mpc/h
        units_override = {
                          'length_unit': ds0.length_unit * ds0.scale_factor,
                          'time_unit'  : ds0.time_unit   * ds0.scale_factor**2,
                          'mass_unit'  : ds0.mass_unit,
                         }

        ds = yt.load( ds0.filename, units_override=units_override )
        del ds0
        gc.collect()

        ds.dataset_index = int(ds.basename[-6:])

        def _particle_index( field, data ):
           return data[('all', 'ParPUID')]
        ds.add_field( ('all', 'particle_index'), function=_particle_index, sampling_type='particle', units='dimensionless' )

        return ds

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]



def check_comoving_units(ds):
    from yt.utilities.cosmology import Cosmology
    co = Cosmology(hubble_constant=ds.parameters['Hubble0'], omega_matter=ds.parameters['OmegaM0'], omega_lambda=1.0-ds.parameters['OmegaM0'])
    a  = ds.scale_factor
    z  = ds.current_redshift

    t_age      = co.t_from_a( a ).in_units('Myr')
    t_lookback = co.lookback_time(0.0, z).in_units('Myr')
    H_of_z     = co.hubble_parameter(z).in_units('km/s/Mpc')

    print( "---------------------------------------------------------" )
    print( f"{ ds.basename                                        = }" )
    print( f"{ ds.current_redshift                                = }" )
    print( f"{ ds.scale_factor                                    = }" )
    print( "" )
    print( f"{ ds.parameters['Hubble0']                           = }" )
    print( f"{ ds.parameters['OmegaM0']                           = }" )
    print( f"{ t_age                                              = }" )
    print( f"{ t_lookback                                         = }" )
    print( f"{ H_of_z                                             = }" )
    print( '' )
    print( f"{ ds.time_unit                                       = }" )
    print( f"{ ds.length_unit                                     = }" )
    print( f"{ ds.mass_unit                                       = }" )
    print( "" )
    print( f"{ ds.time_unit.to('Gyr')                             = }" )
    print( f"{ ds.length_unit.to('Mpc/h')                         = }" )
    print( f"{ ds.mass_unit.to('Msun')                            = }" )
    print( "" )
    print( f"{ ds.domain_width.in_units('code_length')            = }" )
    print( f"{ ds.domain_width.in_units('Mpccm/h')                = }" )
    print( f"{ ds.domain_width.in_units('Mpc/h')                  = }" )
    print( f"{ ds.domain_width.in_units('Mpc')                    = }" )
    print( "" )
    print( f"{ ds.quan(1.0, 'code_length').in_units('Mpccm/h')    = }" )
    print( f"{ ds.quan(1.0, 'code_length').in_units('Mpc/h')      = }" )
    print( f"{ ds.quan(1.0, 'code_time').in_units('Gyr')          = }" )
    print( f"{ ds.quan(1.0, 'code_mass').in_units('Msun')         = }" )
    print( f"{ ds.quan(1.0, 'code_velocity').in_units('km/s')     = }" )
    print( f"{ ds.quan(1.0, 'code_density').in_units('g/cm**3')   = }" )
    print( f"{ ds.quan(1.0, 'code_density').in_units('g/cmcm**3') = }" )
    print( "---------------------------------------------------------" )
    print( "" )

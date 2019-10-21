import yt

# target file
file = 'Data_000000'

# load data
ds = yt.load( file )
ad = ds.all_data()

# calculate center-of-mass
cm = ad.quantities.weighted_average_quantity( ['x','y','z'], weight='cell_mass' )

print( 'CM = (%13.7e, %13.7e, %13.7e) in code units' % (cm[0].in_units('code_length'),
                                                        cm[1].in_units('code_length'),
                                                        cm[2].in_units('code_length')) )

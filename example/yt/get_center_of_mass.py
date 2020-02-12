import yt

# target file
file = 'Data_000000'

# load data
ds = yt.load( file )
ad = ds.all_data()

# filter data if required
# ref: https://yt-project.org/docs/dev/analyzing/filtering.html#cut-regions
dense = ad.cut_region( ['obj["Dens"]>1.0e2'] )

# calculate center-of-mass
cm = dense.quantities.weighted_average_quantity( ['x','y','z'], weight='cell_mass' )

print( 'CM = (%13.7e, %13.7e, %13.7e) in code units' % (cm[0].in_units('code_length'),
                                                        cm[1].in_units('code_length'),
                                                        cm[2].in_units('code_length')) )

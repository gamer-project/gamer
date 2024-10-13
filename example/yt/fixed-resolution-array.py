import yt

# ref: https://yt-project.org/doc/examining/low_level_inspection.html#examining-grid-data-in-a-fixed-resolution-array

# target file
file = 'Data_000000'

# target amr level
lv = 2

# load data
ds = yt.load( file )

# get unsmoothed/smoothed fixed resolution array
ad = ds.covering_grid( level=lv, left_edge=[0.0,0.0,0.0], dims=ds.domain_dimensions*2**lv )
#ad = ds.smoothed_covering_grid( level=lv, left_edge=[0.0,0.0,0.0], dims=ds.domain_dimensions*2**lv )

density = ad['Dens']
print( density.shape )
print( density )

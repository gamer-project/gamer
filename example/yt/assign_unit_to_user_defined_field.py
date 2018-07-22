import yt

ds = yt.load( 'Data_000000' )

# Any user-defined field, for example, "Electron", is by default treated as a
# **dimensionless** YTQuantity. So first define the following function to assign
# the correct code unit to this new field.
# Caveat: do NOT use the same label as the on-disc field (for example, here we use
# "electron" instead of "Electron")

def _electron( field, data ):
   return data["Electron"]*ds.mass_unit/ds.length_unit**3


# Then add this new dimensional field
# Caveat: Here "units" is the output unit, which should be consistent with but
# not necessary the same as the returned unit of _electron()
ds.add_field( ("gamer", "electron"), function=_electron, sampling_type="cell", units="g/cm**3" )

sz = yt.SlicePlot( ds, 'z', 'electron', center='c' )
sz.save()


import numpy as np
import matplotlib.pyplot as plt

Table_Time, Table_PeakError, Table_AveError, Table_MaxError = np.loadtxt( 'DensityProfileErrorEvolution', unpack=True )

# create the figure
fig = plt.figure()
ax  = fig.add_subplot(111)

# plot the profiles
ax.plot( Table_Time, Table_PeakError, label='Center' )
ax.plot( Table_Time, Table_AveError,  label='Ave'    )
ax.plot( Table_Time, Table_MaxError,  label='Max'    )

# set scales
ax.set_xscale('linear')
ax.set_yscale('linear')

# set lables
ax.legend()
ax.set_xlabel( 't'     )
ax.set_ylabel( 'Error' )

ax.set_ylim( 0.0, 1.0 )

# save the figure
plt.tight_layout()
fig.savefig( 'fig_DensityProfile_ErrorEvolution.png' )

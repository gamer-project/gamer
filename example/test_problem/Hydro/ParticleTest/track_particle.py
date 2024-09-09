import yt
import numpy as np
import matplotlib.pylab as plt
import argparse
import sys



#====================================================================================================
# Class
#====================================================================================================
class particle():
    def __init__( self, UID ):
        self.UID    = UID
        self.idx    = 0
        self.t      = []
        self.x_UID  = []
        self.y_UID  = []
        self.z_UID  = []
        self.vx_UID = []
        self.vy_UID = []
        self.vz_UID = []
        self.x_idx  = []
        self.y_idx  = []
        self.z_idx  = []
        self.vx_idx = []
        self.vy_idx = []
        self.vz_idx = []

    def set_idx( self, ds, hdf5=True ):
        if hdf5:
            dd = ds.all_data()
            par_uid  = dd['ParPUid']
            self.idx = np.where(par_uid == self.UID)[0][0]
        else:
            par_uid  = ds[:, 12]
            print(np.where(par_uid == self.UID))
            self.idx = np.where(par_uid == self.UID)[0][0]

    def append_particle_attributes( self, ds, hdf5=True ):
        if hdf5:
            dd = ds.all_data()
            time     = ds.current_time
            posx     = dd['ParPosX']
            posy     = dd['ParPosY']
            posz     = dd['ParPosZ']
            velx     = dd['ParVelX']
            vely     = dd['ParVelY']
            velz     = dd['ParVelZ']
            par_type = dd['ParType']
            par_uid  = dd['ParPUid']
            self.t.append(time)
        else:
            time     = ds[0, 10]
            posx     = ds[:, 1]
            posy     = ds[:, 2]
            posz     = ds[:, 3]
            velx     = ds[:, 4]
            vely     = ds[:, 5]
            velz     = ds[:, 6]
            par_type = ds[:, 11]
            par_uid  = ds[:, 12]
            self.t.append(time)

        self.x_UID.append(posx[par_uid == self.UID])
        self.y_UID.append(posy[par_uid == self.UID])
        self.z_UID.append(posz[par_uid == self.UID])
        self.vx_UID.append(velx[par_uid == self.UID])
        self.vy_UID.append(vely[par_uid == self.UID])
        self.vz_UID.append(velz[par_uid == self.UID])
        self.x_idx.append(posx[self.idx])
        self.y_idx.append(posy[self.idx])
        self.z_idx.append(posz[self.idx])
        self.vx_idx.append(velx[self.idx])
        self.vy_idx.append(vely[self.idx])
        self.vz_idx.append(velz[self.idx])

    def plot_particle( self ):
        lw_UID = 10
        lw_idx = 5
        fig, ax = plt.subplots( 2, 3, sharex='col', figsize=(15, 10) )
        plt.subplots_adjust( hspace=0.1 )

        ax[0, 0].plot( self.t, self.x_UID,  label="UID %05d"%(self.UID), linewidth=lw_UID )
        ax[0, 0].plot( self.t, self.x_idx,  label="idx %05d"%(self.idx), linewidth=lw_idx )
        ax[0, 0].set( ylim=[0, 1] )
        ax[0, 0].tick_params( which="both", direction="in" )

        ax[0, 1].plot( self.t, self.y_UID,  label="UID %05d"%(self.UID), linewidth=lw_UID )
        ax[0, 1].plot( self.t, self.y_idx,  label="idx %05d"%(self.idx), linewidth=lw_idx )
        ax[0, 1].set( ylim=[0, 1] )
        ax[0, 1].tick_params( which="both", direction="in" )

        ax[0, 2].plot( self.t, self.z_UID,  label="UID %05d"%(self.UID), linewidth=lw_UID )
        ax[0, 2].plot( self.t, self.z_idx,  label="idx %05d"%(self.idx), linewidth=lw_idx )
        ax[0, 2].set( ylim=[0, 1] )
        ax[0, 2].tick_params( which="both", direction="in" )

        ax[1, 0].plot( self.t, self.vx_UID, label="UID %05d"%(self.UID), linewidth=lw_UID )
        ax[1, 0].plot( self.t, self.vx_idx, label="idx %05d"%(self.idx), linewidth=lw_idx )
        ax[1, 0].tick_params( which="both", direction="in" )

        ax[1, 1].plot( self.t, self.vy_UID, label="UID %05d"%(self.UID), linewidth=lw_UID )
        ax[1, 1].plot( self.t, self.vy_idx, label="idx %05d"%(self.idx), linewidth=lw_idx )
        ax[1, 1].tick_params( which="both", direction="in" )

        ax[1, 2].plot( self.t, self.vz_UID, label="UID %05d"%(self.UID), linewidth=lw_UID )
        ax[1, 2].plot( self.t, self.vz_idx, label="idx %05d"%(self.idx), linewidth=lw_idx )
        ax[1, 2].tick_params( which="both", direction="in" )

        ax[1, 2].legend()

        ax[0, 0].set(title="x")
        ax[0, 1].set(title="y")
        ax[0, 2].set(title="z")

        ax[1, 0].set(xlabel="time")
        ax[1, 1].set(xlabel="time")
        ax[1, 2].set(xlabel="time")

        ax[0, 0].set(ylabel="pos")
        ax[1, 0].set(ylabel="vel")

        plt.suptitle( "Particle %05d"%self.UID )
        plt.savefig( "particle_%05d.png"%self.UID )
        # plt.show()
        plt.close()
        return


#====================================================================================================
# Main
#====================================================================================================
# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the particles position and velocity by UID and array index' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-u', action='store', required=True, dest='uids', type=int, nargs='+',
                     help='UID of particles' )
parser.add_argument( '-t', action='store_true', required=False, dest='use_hdf5',
                     help='Use Data_* to analysis (slower)' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='data path prefix [%(default)s]', default='./' )

args = parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )

idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix
check_UID = args.uids
use_hdf5  = args.use_hdf5

pars = [ particle(UID) for UID in check_UID ]
for i in range(idx_start, idx_end+1, didx):
    if use_hdf5:
        ds = yt.load( "%s/Data_%06d"%(prefix, i) )
    else:
        print("Loading: %06d"%i)
        ds = np.loadtxt( "%s/Particle_%06d.txt"%(prefix, i) )

    for p in pars:
        if i == 0: p.set_idx( ds, use_hdf5 )

        p.append_particle_attributes( ds, use_hdf5 )

for p in pars:
    p.plot_particle()

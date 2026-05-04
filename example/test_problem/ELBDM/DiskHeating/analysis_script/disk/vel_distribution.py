import matplotlib
matplotlib.use('Agg')
import h5py
import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt


# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
Div_disk = 249        # total number of evenly divided radial bins


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Get the velocity distribution of each 2-kpc bin' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )


#-----------------------------------------------------------------------------------------------------------------
# predefined functions
def SearchIndex(x, A, N):
   i = 0
   j = N - 1
   while(i <= j):
      mid = int(i + (j - i)/2)
      if(A[mid] == x):
         i = mid
         break
      elif(A[mid] > x):
         j = mid - 1
      else: i = mid + 1
   return i

def Gauss(x, mean, disp):
   return np.exp(-(x-mean)**2./2./disp**2.)/disp/(2.*np.pi)**0.5


#-----------------------------------------------------------------------------------------------------------------
# load simulation units and parameters
'''
unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :      100000}
'''
f = h5py.File('../../Data_%06d'%idx_start, 'r')
Unit_L = f['Info']['InputPara']['Unit_L']
Unit_T = f['Info']['InputPara']['Unit_T']*(3.16887646e-14)
Unit_V = f['Info']['InputPara']['Unit_V']
Unit_M = f['Info']['InputPara']['Unit_M']
Center = np.loadtxt('../../Record__Center', skiprows=1, dtype=float)
Center = Center * Unit_L
f.close()
if Center.ndim == 1:
   Center = Center.reshape(1,len(Center)) # reshape the array if there is only one row


#-----------------------------------------------------------------------------------------------------------------
# analyze and plot input simulation data
label = ['$R$ = 4 kpc','$R$ = 6 kpc','$R$ = 8 kpc','$R$ = 10 kpc']
for idx in range(idx_start, idx_end+1, didx):
   print('loading file Data_%06d ...'%idx)
   with h5py.File('../../Data_%06d'%idx, 'r') as f:
      disk_mass = np.array(f['Particle/ParMass'], dtype = np.float64) * Unit_M
      disk_posx = np.array(f['Particle/ParPosX'], dtype = np.float64) * Unit_L
      disk_posy = np.array(f['Particle/ParPosY'], dtype = np.float64) * Unit_L
      disk_posz = np.array(f['Particle/ParPosZ'], dtype = np.float64) * Unit_L
      disk_velx = np.array(f['Particle/ParVelX'], dtype = np.float64) * Unit_V
      disk_vely = np.array(f['Particle/ParVelY'], dtype = np.float64) * Unit_V
      disk_velz = np.array(f['Particle/ParVelZ'], dtype = np.float64) * Unit_V
      disk_type = np.array(f['Particle/ParType'], dtype = np.int32)
      current_step = f['Info/KeyInfo']['Step']
      time = f['Info/KeyInfo']['Time'][0]*Unit_T
   # particle filter: disk particles have ParType=2
   disk_index = (disk_type==2)
   disk_mass = disk_mass[disk_index]
   disk_posx = disk_posx[disk_index]
   disk_posy = disk_posy[disk_index]
   disk_posz = disk_posz[disk_index]
   disk_velx = disk_velx[disk_index]
   disk_vely = disk_vely[disk_index]
   disk_velz = disk_velz[disk_index]
   disk_size = np.size(disk_mass)
   VCM = [ np.sum(disk_mass*disk_velx)/ np.sum(disk_mass),
           np.sum(disk_mass*disk_vely)/ np.sum(disk_mass),
           np.sum(disk_mass*disk_velz)/ np.sum(disk_mass) ]
   center = Center[current_step,3:6]
   disk_posx = disk_posx - center[0]
   disk_posy = disk_posy - center[1]
   disk_posz = disk_posz - center[2]
   disk_velx = disk_velx - VCM[0]
   disk_vely = disk_vely - VCM[1]
   disk_velz = disk_velz - VCM[2]
   # compute angular momentum
   disk_pos = np.zeros((disk_size, 3))
   disk_vel = np.zeros((disk_size, 3))
   disk_pos[:,0] = disk_posx
   disk_pos[:,1] = disk_posy
   disk_pos[:,2] = disk_posz
   disk_vel[:,0] = disk_velx
   disk_vel[:,1] = disk_vely
   disk_vel[:,2] = disk_velz
   disk_L = np.cross(disk_pos, disk_vel)
   tot_L = np.array([np.sum(disk_L[:,0]),np.sum(disk_L[:,1]),np.sum(disk_L[:,2])])
   tot_L = tot_L/(tot_L[0]**2+tot_L[1]**2+tot_L[2]**2)**0.5
   vec = np.cross(tot_L, np.array([0,0,1]))
   c = tot_L[2]
   s = 1-c**2
   matric = np.array([ [0, -vec[2], vec[1]],
                       [vec[2], 0, -vec[0]],
                       [-vec[1], vec[0], 0] ])
   R = np.identity(3) + matric + np.matmul(matric, matric)/(1+c)
   disk_pos_new = np.matmul(R, np.transpose(disk_pos))
   disk_vel_new = np.matmul(R, np.transpose(disk_vel))
   disk_posx = disk_pos_new[0,:]
   disk_posy = disk_pos_new[1,:]
   disk_posz = disk_pos_new[2,:]
   disk_velx = disk_vel_new[0,:]
   disk_vely = disk_vel_new[1,:]
   disk_velz = disk_vel_new[2,:]

   disk_r = (disk_posx**2 + disk_posy**2)**0.5
   disk_sortR  = np.sort(disk_r)/3.08568e+21
   disk_indexR = np.argsort(disk_r)

   num_3 = SearchIndex( 3, disk_sortR, disk_size )
   num_5 = SearchIndex( 5, disk_sortR, disk_size )
   num_7 = SearchIndex( 7, disk_sortR, disk_size )
   num_9 = SearchIndex( 9, disk_sortR, disk_size )
   num_11 = SearchIndex( 11, disk_sortR, disk_size )

   vel = []
   vel.append( np.sort( disk_velz[disk_indexR[num_3:num_5]] )/1e+5 )
   vel.append( np.sort( disk_velz[disk_indexR[num_5:num_7]] )/1e+5 )
   vel.append( np.sort( disk_velz[disk_indexR[num_7:num_9]] )/1e+5 )
   vel.append( np.sort( disk_velz[disk_indexR[num_9:num_11]] )/1e+5 )

   vmean = np.zeros(4)
   disp = np.zeros(4)
   mean_new = np.zeros(4)
   disp_new = np.zeros(4)
   vmax = np.zeros(4)
   vmin = np.zeros(4)
   dv = np.zeros(4)
   index_max = np.zeros(4, dtype = np.int32)
   index_min = np.zeros(4, dtype = np.int32)

   for i in range(4):
      vmean[i] = np.mean(vel[i])
      disp[i] = np.mean( (vel[i] - vmean[i])**2. )**0.5
      vmax[i] = vmean[i] + 3.*disp[i]
      vmin[i] = vmean[i] - 3.*disp[i]
      index_max[i] = SearchIndex( vmax[i], vel[i], len(vel[i]) )
      index_min[i] = SearchIndex( vmin[i], vel[i], len(vel[i]) )
      mean_new[i] = np.mean(vel[i][index_min[i]:index_max[i]])
      disp_new[i] = np.mean( (vel[i][index_min[i]:index_max[i]] - mean_new[i])**2. )**0.5
      dv[i] = (vmax[i]-vmin[i])/Div_disk

   Data = np.zeros((4, 2, Div_disk))
   for i in range(4):
      num_v = index_min[i]
      v = vmin[i]
      for j in range(Div_disk):
         num_v_pre = num_v
         Data[i,0,j] = v + 0.5*dv[i]
         num_v = SearchIndex( v + dv[i], vel[i], len(vel[i]) )
         Data[i,1,j] = (num_v - num_v_pre)/dv[i]/(index_max[i]- index_min[i])
         v = v + dv[i]

   np.savez('Vel_data_%06d'%idx, data = Data, mean = mean_new, disp = disp_new, time_Myr = time)

   print('plotting vel_%06d.png ...'%idx)
   fig, axs = plt.subplots(2, 2)
   fig_num = 0
   for i in range(2):
      for j in range(2):
         axs[i, j].plot(Data[fig_num,0,:], Data[fig_num,1,:], color = 'blue')
         axs[i, j].plot(Data[fig_num,0,:], Gauss(Data[fig_num,0,:],mean_new[fig_num],disp_new[fig_num]),'--', lw=0.8, color = 'red', label ='Gaussian')
         axs[i, j].set_title(label[fig_num], fontsize=6)
         axs[i, j].set_xlim((vmean[fig_num]-3*disp[fig_num],vmean[fig_num]+3*disp[fig_num]))
         axs[i, j].tick_params(axis='both', labelsize=6)
         axs[i, j].grid(ls='--',lw=0.3)
         axs[i, j].legend(loc='best', shadow=True, fontsize=6)

         fig_num += 1
   axs[0, 0].set_ylabel('number density', fontsize=6)
   axs[1, 0].set_ylabel('number density', fontsize=6)
   axs[1, 0].set_xlabel('$v(km/s)$', fontsize=6 )
   axs[1, 1].set_xlabel('$v(km/s)$', fontsize=6 )

   fig.tight_layout()
   fig.suptitle('t = %4.1f Myr' %time, fontsize=8)
   plt.savefig("vel_%06d.png"%idx, dpi = 140)
   plt.close()

   print(' ')

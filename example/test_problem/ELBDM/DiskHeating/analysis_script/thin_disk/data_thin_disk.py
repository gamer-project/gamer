import h5py
import math
import numpy as np
import argparse
import sys


# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
Div_disk = 500             # total number of data points sampled per disk rotation curve
ratio    = 0.76159415595   # tanh(1)


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Get disk properties' )

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
# analyze simulation data and generate processed output files
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
   # particle filter: thin disk particles have ParType=3
   disk_index = (disk_type==3)
   disk_mass = disk_mass[disk_index]
   disk_posx = disk_posx[disk_index]
   disk_posy = disk_posy[disk_index]
   disk_posz = disk_posz[disk_index]
   disk_velx = disk_velx[disk_index]
   disk_vely = disk_vely[disk_index]
   disk_velz = disk_velz[disk_index]
   disk_size = np.size(disk_mass)
   print("number of thin disk particles = %d"%disk_size)
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
   disk_velr = (disk_posx*disk_velx + disk_posy*disk_vely)/disk_r
   disk_velp = (disk_posx*disk_vely - disk_posy*disk_velx)/disk_r
   disk_sortR  = np.sort(disk_r)
   disk_indexR = np.argsort(disk_r)
   Data = np.zeros((8,Div_disk))

   r = 0
   disk_num = 0
   dr = 15.0*3.08568e+21/Div_disk
   for j in range(Div_disk):
      disk_num_pre = disk_num
      Data[0,j] = r + 0.5 * dr
      disk_num = SearchIndex( r+dr, disk_sortR, disk_size )
      mean_vr = np.mean( disk_velr[ disk_indexR[ disk_num_pre:disk_num ] ] )
      mean_vp = np.mean( disk_velp[ disk_indexR[ disk_num_pre:disk_num ] ] )
      mean_vz = np.mean( disk_velz[ disk_indexR[ disk_num_pre:disk_num ] ] )
      mean_vp2 = np.mean( disk_velp[ disk_indexR[ disk_num_pre:disk_num ] ]**2 )
      Data[1,j] = mean_vp      # rotation curve
      Data[2,j] = (np.mean(( disk_velr[ disk_indexR[ disk_num_pre:disk_num ] ] - mean_vr )**2))**0.5      # sigma_r
      Data[3,j] = (np.mean(( disk_velz[ disk_indexR[ disk_num_pre:disk_num ] ] - mean_vz )**2))**0.5      # sigma_z
      mass = disk_mass[0]*( disk_num-disk_num_pre )
      mass_total = disk_mass[0]*disk_num
      area = (math.pi *((r + dr)**2 - r**2))
      Data[4,j] = mass/area    # surface density Sigma
      Data[5,j] = mean_vp2     # sigma_phi^2
      Data[6,j] = mass_total   # enclosed mass
      # get sacle height
      CMZ = np.average(disk_posz[ disk_indexR[ disk_num_pre:disk_num ] ])
      disk_z = np.abs( disk_posz[ disk_indexR[ disk_num_pre:disk_num ] ] - CMZ )
      sortZ  = np.sort(disk_z)
      target_index = ratio*len(disk_z)
      index_plus = int(target_index)
      index_minus = int(target_index - 1)
      target_z = sortZ[index_minus] + (target_index - index_minus)*(sortZ[index_plus] - sortZ[index_minus])
      Data[7,j] = target_z     # scale height

      r = r + dr

   np.save('Data_Thin_Disk_%06d'%idx, Data)
   print(' ')

import h5py
import numpy as np

filename = "snap_001.hdf5"
UNIT_L_GALIC = 3.08568e+21
UNIT_M_GALIC = 1.989e+43
UNIT_D_GALIC = UNIT_M_GALIC/UNIT_L_GALIC**3
UNIT_V_GALIC = 1.0e+5

UNIT_L = 4.436632034507548e+24
UNIT_D = 2.5758579703106658e-30
UNIT_M = UNIT_D*UNIT_L**3
UNIT_V = 1.0e+7

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

with h5py.File(filename, "r") as f:
    N_Disk = f['PartType2']['Masses'].size
    mass = f['PartType2']['Masses'][0] *UNIT_M_GALIC/UNIT_M
    disk_posx = np.array(f['PartType2']['Coordinates'][:,0], dtype = np.float64) *UNIT_L_GALIC/UNIT_L
    disk_posy = np.array(f['PartType2']['Coordinates'][:,1], dtype = np.float64) *UNIT_L_GALIC/UNIT_L
    disk_posz = np.array(f['PartType2']['Coordinates'][:,2], dtype = np.float64) *UNIT_L_GALIC/UNIT_L
    disk_velx = np.array(f['PartType2']['Velocities'][:,0], dtype = np.float64) *UNIT_V_GALIC/UNIT_V
    disk_vely = np.array(f['PartType2']['Velocities'][:,1], dtype = np.float64) *UNIT_V_GALIC/UNIT_V
    disk_velz = np.array(f['PartType2']['Velocities'][:,2], dtype = np.float64) *UNIT_V_GALIC/UNIT_V

CM = [ np.sum(disk_posx)/ N_Disk,
       np.sum(disk_posy)/ N_Disk,
       np.sum(disk_posz)/ N_Disk ]
VCM = [ np.sum(disk_velx)/ N_Disk,
        np.sum(disk_vely)/ N_Disk,
        np.sum(disk_velz)/ N_Disk ]
#print(np.array(VCM)*UNIT_V/UNIT_V_GALIC)

center = [7.9770907e-02, 7.9854590e-02, 8.0858787e-02]
disk_posx = disk_posx - CM[0]
disk_posy = disk_posy - CM[1]
disk_posz = disk_posz - CM[2]
disk_velx = disk_velx - VCM[0]
disk_vely = disk_vely - VCM[1]
disk_velz = disk_velz - VCM[2]
disk_r = (disk_posx**2+disk_posy**2+disk_posz**2)**0.5
sortR_kpc = np.sort(disk_r)*UNIT_L/UNIT_L_GALIC

#print(sortR_kpc[0])
#print(sortR_kpc[-1])

indexR = np.argsort(disk_r)

# exclude particles with r > 70 kpc
num = SearchIndex( 70, sortR_kpc, N_Disk )

disk_posx = disk_posx[ indexR[:num] ] + center[0]
disk_posy = disk_posy[ indexR[:num] ] + center[1]
disk_posz = disk_posz[ indexR[:num] ] + center[2]
disk_velx = disk_velx[ indexR[:num] ]
disk_vely = disk_vely[ indexR[:num] ]
disk_velz = disk_velz[ indexR[:num] ]

disk_type = np.full(num, 2)
disk_mass = np.full(num, mass)
with open('DiskHeatingParticleIC', 'wb') as output:
    output.write(disk_mass.astype('f').tobytes())
    output.write(disk_posx.astype('f').tobytes())
    output.write(disk_posy.astype('f').tobytes())
    output.write(disk_posz.astype('f').tobytes())
    output.write(disk_velx.astype('f').tobytes())
    output.write(disk_vely.astype('f').tobytes())
    output.write(disk_velz.astype('f').tobytes())
    output.write(disk_type.astype('f').tobytes())
output.close()
#print('Disk particles = '+ str(len(disk_posx)))
print('Disk particles = '+ str(num))
print('DiskHeatingParticleIC complete')

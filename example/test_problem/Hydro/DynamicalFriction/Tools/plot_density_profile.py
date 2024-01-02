import argparse
import sys
import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

center_mode = 'max'
dpi         = 300
nbin        = 128

def NFW_func(r, rho_0, R_s):
     return rho_0/(r/R_s*(1+r/R_s)**2)
NFW = [1.96311379e+06,5.48808202e+00]
x = np.linspace(0.001,100,1000)

def Burkert_dens(x):
    return ((1/((1+x)*(1+x*x))))
def Plummer_dens(x):
    return (1+x*x)**(-2.5)
def density(rho_s,r0,r,model_name):
    x  =  r/r0
    if model_name =="Burkert":
        return rho_s*Burkert_dens(x)
    elif model_name =="Plummer":
        return rho_s*Plummer_dens(x)
    else:
        return 0
y_Burkert = density(562.63979048052,0.25,x,"Burkert")*1E6
y_Plummer = density(37.12790416639654,0.668,x,"Plummer")*1E6
y_sum = y_Burkert+y_Plummer

yt.enable_parallelism()
#ts = yt.load( [ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
files = [ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ]
count = idx_start

ds0 = yt.load("../Data_000000")


for f in files:
   
   ds = yt.load(f)
   time = round(ds.current_time,3)
   def find_GC(pfilter,data):
        filter = data[('all', u'ParType')] == 4
        return filter
   yt.add_particle_filter("GC",function=find_GC,filtered_type="all",requires=["ParType"])

   ds.add_particle_filter("GC")
   ad = ds.all_data()
   min_potential, min_position = ds.find_min(('gamer',u'Pote'))
   center_halo = min_position.value 
   sp1 = ds.sphere(center_halo, (10000, "code_length")) #ds.sphere( center_mode, 0.5*ds.domain_width.to_value().max() )
   sp0 = ds0.sphere(center_halo,(10000,"code_length"))

   field ='particle_density_on_grid'
   
   rp1 = yt.create_profile(sp1,"radius",field,units={("radius"): "kpc"},logs={("radius"): True},n_bins=nbin,extrema = {'radius': ((0.01, 'kpc'), (15, 'kpc'))})
   rp0 = yt.create_profile(sp0,"radius",field,units={("radius"): "kpc"},logs={("radius"): True},n_bins=nbin,extrema = {'radius': ((0.01, 'kpc'), (15, 'kpc'))})
   
   fig = plt.figure(dpi=dpi)
   print('plotting....')
   ax = fig.add_subplot(111)
   
   ax.plot(rp0.x.value,rp0[field].in_units("msun/kpc**3").value,label='DATA_000000',color='blue',linewidth=0.75)
   ax.plot(rp1.x.value,rp1[field].in_units("msun/kpc**3").value,label='',color='blue',linewidth=0.75)
#   ax.plot(x,NFW_func(x,NFW[0],NFW[1]),label='NFW log fit',color='#ff3c00',linewidth=0.5)

   
#   df = pd.read_csv('../profile_Fornax.txt',sep=' ',skiprows=1,header=None)
#   ax.plot(df[0],df[1]*1E6,label='profile_Fornax.txt',alpha=0.5,color='red',linewidth=0.75)
#   ax.set_ylim(10,5E7)
#   ax.plot(x,y_Burkert,label='Halo Part',color='#DCBACE',linewidth=0.5)
#   ax.plot(x,y_Plummer,label='Stellar Part',color='#BACEDC',linewidth=0.5)
#   ax.plot(x,y_sum,label='Halo+Stellar',color='#cedcba',linewidth=0.5)

   ax.set_xlim(0.01,10)
   ax.set_xscale('log')
   ax.set_yscale('log')
   fig.text(0.5,0.04,'Radius(kpc)',ha='center')
   fig.text(0.02,0.5,r'Density($\dfrac{M_{\odot}}{kpc^3}$)',va='center',rotation='vertical')
   ax.legend(loc='best')
   #ax.set_text("Burkert fit N_PAR=100000")
   fig.savefig(f+"_profile.png")
   print(f+" DONE!!!!!")
   count+=1

   

#for f in files:
#
#  ds = yt.load(f)
#   time = round(ds.current_time * 14.06846789227,3)
#   sp1 = ds.sphere("c", (1.5, "code_length")) #ds.sphere( center_mode, 0.5*ds.domain_width.to_value().max() )
#   sp0 = d0.sphere("c", (1.5, "code_length"))
#   field ='particle_density_on_grid'
#   rp0 = yt.create_profile(sp0,"radius",field,units={("radius"): "kpc"},logs={("radius"): True},n_bins=128,extrema = {'radius': ((0.05, 'kpc'), (105, 'kpc'))})
#   rp1 = yt.create_profile(sp1,"radius",field,units={("radius"): "kpc"},logs={("radius"): True},n_bins=128,extrema = {'radius': ((0.05, 'kpc'), (105, 'kpc'))})
#   fig = plt.figure(dpi=dpi)
#   print('plotting....')
#   ax = fig.add_subplot(111)
#   ax.plot(rp0.x.value,rp0[field].in_units("msun/kpc**3").value,label='0 Gyr',linewidth=0.8)
#   ax.plot(rp1.x.value,rp1[field].in_units("msun/kpc**3").value,label=str(time)+' Gyr',color='blue',linewidth=0.8)
#   ax.plot(x,NFW_func(x,NFW[0],NFW[1]),label='NFW log fit',color='#ff3c00',linewidth=0.5)
#   ax.set_ylim(10,5E7)
#   ax.set_xlim(0.5,100)
#   ax.set_xscale('log')
#   ax.set_yscale('log')
#   fig.text(0.5,0.04,'Radius(kpc)',ha='center')
#   fig.text(0.02,0.5,r'Density($\dfrac{M_{\odot}}{kpc^3}$)',va='center',rotation='vertical')
#   ax.legend()
#   #ax.set_text("NFW fit N_PAR=100000 ")
#   fig.savefig(f+"_profile.png")
#   print(f+" Done!!!")

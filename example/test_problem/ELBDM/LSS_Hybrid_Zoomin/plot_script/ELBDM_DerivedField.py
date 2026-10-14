import yt
import numpy as np

# define the derived fields
###########################
# psi = R + iI = f*e^(iS)

# amplitude
# = sqrt(rho) = sqrt(R^2 + I^2)
def _f( field, data ):
   return data["Dens"]**0.5

# real part
# = sqrt(rho)*cos(S) = R
def _Real( field, data ):
   return data["f"]*np.cos( data["Phase"] )

# imaginary part
# = sqrt(rho)*sin(S) = I
def _Imag( field, data ):
   return data["f"]*np.sin( data["Phase"] )

# phase
# = arctan(I/R)
def _S( field, data ):
   return np.arctan2( data["Imag"], data["Real"] )

def _phasemod2( field, data ):
   return np.arctan2( np.sin(data["Phase"]), np.cos(data["Phase"]) )

# momentum density x-component
# = rho*v_x = (R*dI/dx - I*dR/dx)*hbar/m
def _j_x( field, data ):
   ELBDM_ETA = data.ds.parameters['ELBDM_Mass']*data.ds.units.code_mass/data.ds.units.reduced_planck_constant
   return (data["Real"]*data["Imag_gradient_x"] - data["Imag"]*data["Real_gradient_x"])/ELBDM_ETA

# momentum density y-component
# = rho*v_y = (R*dI/dy - I*dR/dy)*hbar/m
def _j_y( field, data ):
   ELBDM_ETA = data.ds.parameters['ELBDM_Mass']*data.ds.units.code_mass/data.ds.units.reduced_planck_constant
   return (data["Real"]*data["Imag_gradient_y"] - data["Imag"]*data["Real_gradient_y"])/ELBDM_ETA

# momentum density z-component
# = rho*v_z = (R*dI/dz - I*dR/dz)*hbar/m
def _j_z( field, data ):
   ELBDM_ETA = data.ds.parameters['ELBDM_Mass']*data.ds.units.code_mass/data.ds.units.reduced_planck_constant
   return (data["Real"]*data["Imag_gradient_z"] - data["Imag"]*data["Real_gradient_z"])/ELBDM_ETA

# momentum density magnitude
# = rho*|v| = |R*grad(I) - I*grad(R)|*hbar/m
def _j( field, data ):
   return (data["j_x"]**2 + data["j_y"]**2 + data["j_z"]**2)**0.5

# bulk velocity x-component
# = (dS/dx)*hbar/m = (R*dI/dx - I*dR/dx)/rho*hbar/m
def _v_x( field, data ):
   return data["j_x"]/data["Dens"]

# bulk velocity y-component
# = (dS/dy)*hbar/m = (R*dI/dy - I*dR/dy)/rho*hbar/m
def _v_y( field, data ):
   return data["j_y"]/data["Dens"]

# bulk velocity z-component
# = (dS/dz)*hbar/m = (R*dI/dz - I*dR/dz)/rho*hbar/m
def _v_z( field, data ):
   return data["j_z"]/data["Dens"]

# bulk velocity magnitude
# = |grad(S)|*hbar/m = |(R*grad(I) - I*grad(R))|/rho*hbar/m
def _v( field, data ):
   return (data["v_x"]**2 + data["v_y"]**2 + data["v_z"]**2)**0.5

###########################


def Add_derived_fields( ds ):

   ds.add_field( ("gamer","f"),
                 function=_f, units="code_mass**0.5/code_length**1.5",
                 display_name=r"$f$", sampling_type="cell" )
   ds.add_field( ("gamer","Real"),
                 function=_Real, units="code_mass**0.5/code_length**1.5",
                 display_name=r"$Real$", sampling_type="cell" )
   ds.add_field( ("gamer","Imag"),
                 function=_Imag, units="code_mass**0.5/code_length**1.5",
                 display_name=r"$Imag$", sampling_type="cell" )
   ds.add_field( ("gamer","S"),
                 function=_S, units="dimensionless",
                 display_name=r"$S$", sampling_type="cell" )
   ds.add_field( ("gamer", "PhaseMod2"),
                 function=_phasemod2,
                 sampling_type="local", units="" )

   # gradient fields
   Grad_R = ds.add_gradient_fields( ("gamer","Real") )
   Grad_I = ds.add_gradient_fields( ("gamer","Imag") )
   Grad_f = ds.add_gradient_fields( ("gamer","f") )
   Grad_S = ds.add_gradient_fields( ("gamer","S") )

   # momentum density
   ds.add_field( ("gamer","j_x"),
                 function=_j_x, units="code_mass/(code_length**2*code_time)",
                 display_name=r"$j_x$", sampling_type="cell" )
   ds.add_field( ("gamer","j_y"),
                 function=_j_y, units="code_mass/(code_length**2*code_time)",
                 display_name=r"$j_y$", sampling_type="cell" )
   ds.add_field( ("gamer","j_z"),
                 function=_j_z, units="code_mass/(code_length**2*code_time)",
                 display_name=r"$j_z$", sampling_type="cell" )
   ds.add_field( ("gamer","j"),
                 function=_j, units="code_mass/(code_length**2*code_time)",
                 display_name=r"$|j|$", sampling_type="cell" )

   # velocity
   ds.add_field( ("gamer","v_x"),
                 function=_v_x, units="code_length/code_time",
                 display_name=r"$v_x$", sampling_type="cell" )
   ds.add_field( ("gamer","v_y"),
                 function=_v_y, units="code_length/code_time",
                 display_name=r"$v_y$", sampling_type="cell" )
   ds.add_field( ("gamer","v_z"),
                 function=_v_z, units="code_length/code_time",
                 display_name=r"$v_z$", sampling_type="cell" )
   ds.add_field( ("gamer","v"),
                 function=_v, units="code_length/code_time",
                 display_name=r"$|v|$", sampling_type="cell" )

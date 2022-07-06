import yt
import numpy as np
from scipy import constants

# yt reference
# 1. gradient fields
#    https://yt-project.org/doc/analyzing/fields.html?highlight=gradient#gradient-fields 
#    example:
#    https://yt-project.org/doc/cookbook/calculating_information.html#complicated-derived-fields
# 2. weight average
#    https://yt-project.org/doc/reference/api/yt.data_objects.derived_quantities.html?highlight=center%20mass#yt.data_objects.derived_quantities.WeightedAverageQuantity

def momentum_field_x(field, data):
   grad_fields_real_x = data["gamer", "Real_gradient_x"]
   grad_fields_imag_x = data["gamer", "Imag_gradient_x"]
   momentum_field_x   = ( data["gamer","Real"]*grad_fields_imag_x - data["gamer","Imag"]*grad_fields_real_x )/data["gamer","Dens"]
   return momentum_field_x

def momentum_field_y(field, data):
   grad_fields_real_y = data["gamer", "Real_gradient_y"]
   grad_fields_imag_y = data["gamer", "Imag_gradient_y"]
   momentum_field_y   = ( data["gamer","Real"]*grad_fields_imag_y - data["gamer","Imag"]*grad_fields_real_y )/data["gamer","Dens"]
   return momentum_field_y

def momentum_field_z(field, data):
   grad_fields_real_z = data["gamer", "Real_gradient_z"]
   grad_fields_imag_z = data["gamer", "Imag_gradient_z"]
   momentum_field_z   = ( data["gamer","Real"]*grad_fields_imag_z - data["gamer","Imag"]*grad_fields_real_z )/data["gamer","Dens"]
   return momentum_field_z

m_22             = 0.2
h_0              = 0.6732117
r_c              = 5.38569360e-04  # soliton radius, Mpc/h
radius_factor    = 2.0             # use a sphere with radius equals radius_factor*r_c to define soliton
kpc_to_meter     = 3.08567758e19   # kpc -> meter

ds               = yt.load("./RESTART_RESET")
p_x, p_y, p_z    = ("gamer","momentum_field_x"), ("gamer","momentum_field_y"), ("gamer","momentum_field_z")
grad_fields_real = ds.add_gradient_fields(("gamer", "Real"))
grad_fields_imag = ds.add_gradient_fields(("gamer", "Imag"))

ds.add_field(p_x, function=momentum_field_x, sampling_type="cell", units="1/code_length")
ds.add_field(p_y, function=momentum_field_y, sampling_type="cell", units="1/code_length")
ds.add_field(p_z, function=momentum_field_z, sampling_type="cell", units="1/code_length")

source           = ds.sphere( "m", (radius_factor*r_c, "Mpc/h"))
p_x_mean         = np.array(source.quantities.weighted_average_quantity(p_x, weight="cell_mass").to("1/code_length"))  # h/Mpc, computed by \int (real*\nabla_x imag - imag*\nabla_x real)dV/(int \rho*dv)
p_y_mean         = np.array(source.quantities.weighted_average_quantity(p_y, weight="cell_mass").to("1/code_length"))  # h/Mpc, computed by \int (real*\nabla_y imag - imag*\nabla_y real)dV/(int \rho*dv)
p_z_mean         = np.array(source.quantities.weighted_average_quantity(p_z, weight="cell_mass").to("1/code_length"))  # h/Mpc, computed by \int (real*\nabla_z imag - imag*\nabla_z real)dV/(int \rho*dv)

p_x_mean         *= constants.hbar/(m_22*1e-22*constants.e/constants.c**2)*h_0/1e3/kpc_to_meter
p_y_mean         *= constants.hbar/(m_22*1e-22*constants.e/constants.c**2)*h_0/1e3/kpc_to_meter
p_z_mean         *= constants.hbar/(m_22*1e-22*constants.e/constants.c**2)*h_0/1e3/kpc_to_meter

print("Soliton CoM velocity is (%.8e,%.8e,%.8e) m/s"%(p_x_mean, p_y_mean, p_z_mean))

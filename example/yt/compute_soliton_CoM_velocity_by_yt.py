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

### define derived fields for wake vector field k_{x,y,z} with unit specified in add_field
def wave_vector_field_x(field, data):
   grad_fields_real_x    = data["gamer", "Real_gradient_x"]
   grad_fields_imag_x    = data["gamer", "Imag_gradient_x"]
   wave_vector_field_x = ( data["gamer","Real"]*grad_fields_imag_x - data["gamer","Imag"]*grad_fields_real_x )/data["gamer","Dens"]
   return wave_vector_field_x

def wave_vector_field_y(field, data):
   grad_fields_real_y    = data["gamer", "Real_gradient_y"]
   grad_fields_imag_y    = data["gamer", "Imag_gradient_y"]
   wave_vector_field_y = ( data["gamer","Real"]*grad_fields_imag_y - data["gamer","Imag"]*grad_fields_real_y )/data["gamer","Dens"]
   return wave_vector_field_y

def wave_vector_field_z(field, data):
   grad_fields_real_z = data["gamer", "Real_gradient_z"]
   grad_fields_imag_z = data["gamer", "Imag_gradient_z"]
   wave_vector_field_z   = ( data["gamer","Real"]*grad_fields_imag_z - data["gamer","Imag"]*grad_fields_real_z )/data["gamer","Dens"]
   return wave_vector_field_z
###

m_22             = 0.2             # m_a/1e-22
h_0              = 0.6732117       # H_0/100
r_c              = 5.38569360e-04  # soliton radius, Mpc/h
radius_factor    = 2.0             # use a sphere with radius equals radius_factor*r_c to define soliton
kpc_to_meter     = 3.08567758e19   # kpc -> meter

ds               = yt.load("./RESTART_RESET")
k_x, k_y, k_z    = ("gamer","wave_vector_field_x"), ("gamer","wave_vector_field_y"), ("gamer","wave_vector_field_z")
grad_fields_real = ds.add_gradient_fields(("gamer", "Real"))
grad_fields_imag = ds.add_gradient_fields(("gamer", "Imag"))

### add wave vector field k_{x,y,z} in unit of code_length^{-1}, i.e. (Mpc/h)^{-1} here
ds.add_field(k_x, function=wave_vector_field_x, sampling_type="cell", units="1/code_length")
ds.add_field(k_y, function=wave_vector_field_y, sampling_type="cell", units="1/code_length")
ds.add_field(k_z, function=wave_vector_field_z, sampling_type="cell", units="1/code_length")
###

### compute cell_mass weighted wave vector k_{x,y,z} in unit of code_length^{-1}, by \int (real*\nabla_{x,y,z} imag - imag*\nabla_{x,y,z} real)dV/(int \rho*dv)
source           = ds.sphere( "m", (radius_factor*r_c, "Mpc/h"))
k_x_mean         = np.array(source.quantities.weighted_average_quantity(k_x, weight="cell_mass").to("1/code_length"))
k_y_mean         = np.array(source.quantities.weighted_average_quantity(k_y, weight="cell_mass").to("1/code_length"))
k_z_mean         = np.array(source.quantities.weighted_average_quantity(k_z, weight="cell_mass").to("1/code_length"))
###

### compute cell_mass weighted velocity field v_{x,y,z} in unit of m/s, by h_bar*k_{x,y,z}/m_a
v_x_mean         = k_x_mean*constants.hbar/(m_22*1e-22*constants.e/constants.c**2)*h_0/1e3/kpc_to_meter
v_y_mean         = k_y_mean*constants.hbar/(m_22*1e-22*constants.e/constants.c**2)*h_0/1e3/kpc_to_meter
v_z_mean         = k_z_mean*constants.hbar/(m_22*1e-22*constants.e/constants.c**2)*h_0/1e3/kpc_to_meter
###

print("Soliton CoM velocity is (%.8e,%.8e,%.8e) m/s"%(v_x_mean, v_y_mean, v_z_mean))

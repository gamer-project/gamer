"""
This is for analytical solution of shock tube test.
Solve variables:
rho2, rho3, rho4, p2, p3, p4, v2, v3, v4

Shock tube:
--------------------
Left 5|4|3|2|1 Right
--------------------
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

# Fix parameter
gamma_th = 5./3.
gamma_cr = 4./3.
delta_gamma = gamma_th - gamma_cr

# Main function
def shock_sol(rhoL, rhoR, p_crL, p_crR, p_thL, p_thR, vL, vR, x, center, t):
	# define some functions
	# Incomplete beta function
	def B(x, a, b):
		def F(t):
			if t == 0:
				return 0.
			elif t == 1:
				return 1.
			else:
				return t**(a-1) * (1-t)**(b-1)

		sol = quad(F, 0., x)
		return sol[0]

	# I function
	def I(rho, p_th, p_cr):
		x_rho  = gamma_th * p_th / (gamma_cr * p_cr + gamma_th * p_th)
		first  = np.sqrt(gamma_cr * p_cr * rho**(-gamma_cr)) / delta_gamma 
		second = ( gamma_cr * p_cr * rho**(-gamma_cr) / gamma_th / p_th / rho**(-gamma_th) )**((gamma_cr - 1) / (2 * delta_gamma))
		third  = B(x_rho, (gamma_cr - 1)/2/delta_gamma, (1 - gamma_th)/2/delta_gamma)
		return first * second * third

	# Setup the density equations to solve rho2 and rho3
	def sol_rho(z):
		xs = z[0] / rhoR
		xr = z[1] / rhoL
		p_cr2 = p_crR * xs**gamma_cr
		p2    = p_crL * xr**gamma_cr + p_thL * xr**gamma_th
		pR    = p_thR + p_crR
		e_cr2 = p_crR * xs**(gamma_cr) / (gamma_cr - 1)
		e_th2 = (p2 - p_cr2) / (gamma_th - 1)
		p_th3 = p_thL * xr**gamma_th
		p_cr3 = p_crL * xr**gamma_cr

		eq1 = (p2 - pR)*(xs - 1) - rhoR * xs * (I(rhoL, p_thL, p_crL) - I(xr*rhoL, p_th3, p_cr3))**2
		eq2 = (p2 + pR)*(xs - 1) + 2 * (xs * (p_crR/(gamma_cr-1) + p_thR/(gamma_th-1)) - e_th2 - e_cr2)

		return [eq1, eq2]
	
	# Solve the variables
	rho_guess = [4., 0.5]
	temp  = fsolve(sol_rho, rho_guess)
	rho2  = temp[0]
	rho3  = temp[1]

	p_cr2 = p_crR * (rho2/rhoR)**gamma_cr
	p_cr3 = p_crL * (rho3/rhoL)**gamma_cr
	p_th3 = p_thL * (rho3/rhoL)**gamma_th
	p_th2 = p_cr3 + p_th3 - p_cr2
	

	v2    = np.sqrt((p_cr2 + p_th2 - p_crR - p_thR) * (rho2 - rhoR) / rho2 / rhoR) 
	v3    = v2
	vs    = rho2 * v2 / (rho2 - rhoR)
	vt    = I(rho3, p_th3, p_cr3) - I(rhoL, p_thL, p_crL) + np.sqrt(gamma_cr * p_cr3 / rho3 + gamma_th * p_th3 / rho3)
	v5 = np.sqrt((gamma_cr * p_crL + gamma_th * p_thL) / rhoL)

	rho4 = 0.0
	p_cr4 = 0.0
	p_th4 = 0.0
	v4 = 0.0

	# Set the value to the output
	
	def sol_rho4(rho4, x, t):
		#global p_thL, p_crL
		A_th4 = gamma_th * p_thL / rhoL**gamma_th
		A_cr4 = gamma_cr * p_crL / rhoL**gamma_cr
		def I4(rho):
			x_rho = A_th4 * rho**gamma_th / (A_th4 * rho**gamma_th + A_cr4 * rho**gamma_cr)
			first  = np.sqrt(A_cr4) / delta_gamma 
			second = (A_cr4 / A_th4)**((gamma_cr - 1)/2/delta_gamma)
			third  = B(x_rho, (gamma_cr - 1)/2/delta_gamma, (1 - gamma_th)/2/delta_gamma)
			return first * second * third
		eq = I4(rho4) - I(rhoL, p_thL, p_crL) + x/t + np.sqrt(A_cr4 * rho4**(gamma_cr - 1) + A_th4 * rho4**(gamma_th - 1))
		return eq

	L = len(x)
	rho  = np.zeros(L)
	p_cr = np.zeros(L)
	p_th = np.zeros(L)
	vel  = np.zeros(L)

	for i in range(L):
		if x[i] < center - v5 * t:
			rho[i]  = rhoL
			p_cr[i] = p_crL
			p_th[i] = p_thL
			vel[i]  = vL
		elif x[i] < center - vt * t:
			rho[i]  = fsolve(lambda k:sol_rho4(k, x[i] - center, t), 0.5)[0]
			A_th4 = gamma_th * p_thL / rhoL**gamma_th
			A_cr4 = gamma_cr * p_crL / rhoL**gamma_cr
			p_cr[i] = A_cr4 * rho[i]**gamma_cr / gamma_cr
			p_th[i] = A_th4 * rho[i]**gamma_th / gamma_th
			vel[i]  = (x[i] - center)/t + np.sqrt(A_cr4 * rho[i]**(gamma_cr - 1) + A_th4 * rho[i]**(gamma_th - 1))
		elif x[i] < center + v2 * t:
			rho[i]  = rho3
			p_cr[i] = p_cr3
			p_th[i] = p_th3
			vel[i]  = v3
		elif x[i] < center + vs * t:
			rho[i]  = rho2
			p_cr[i] = p_cr2
			p_th[i] = p_th2
			vel[i]  = v2
		else:
			rho[i]  = rhoR
			p_cr[i] = p_crR
			p_th[i] = p_thR
			vel[i]  = vR
	
	return rho, p_cr, p_th, vel

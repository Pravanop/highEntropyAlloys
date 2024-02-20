import matplotlib.pyplot as plt
import numpy as np

def hmix_plotter(x_dft , y_dft , x , y) :
	"""

	:param x_dft:
	:param y_dft:
	:param x:
	:param y:
	"""
	plt.plot(x , y , c = "#882255" , alpha = 0.7 , linewidth = 4)
	plt.scatter(x_dft , y_dft , s = 100 , c = "#882255" , edgecolors = "#000000")
	plt.xlabel("mole fraction")
	plt.ylabel("Mixing Enthalpy (eV/atom)")
	
def phase_diagram_plotter(tangent_points , system) :
	"""

	:param tangent_points:
	:param system:
	"""
	xA , yA = tangent_points[: , 0 , 0] , tangent_points[: , 0 , 1]
	xB , yB = tangent_points[: , 1 , 0] , tangent_points[: , 1 , 1]
	
	x_phase = np.concatenate((xA , np.flip(xB)))
	y = np.concatenate((yA , np.flip(yB)))
	plt.plot(x_phase , y , linewidth = 2 , linestyle = "--" , c = "#882255")
	# plt.text(system[0] , x = 0)
	# plt.text(system[1] , x = 1)
	plt.ylabel('T (K)')
	plt.xlabel('mol_fraction')

def gibbs_vs_T_plotter(gibbs_T , step, T_range) :
	"""

	:param gibbs_T:
	:param step:
	:param T_range:
	"""
	for i in range(0 , len(gibbs_T) + 1 , int(step / 75)) :
		plt.plot(x , gibbs_T[i] , linewidth = 2 , c = "#882255" , alpha = (i + 20) / 100)
	
	plt.legend(range(T_range[0] , T_range[1] , step))
	plt.ylabel("Gibbs Free Energy (meV/atom)")
	plt.xlabel("mol fraction")
	plt.axhline(y = 0)
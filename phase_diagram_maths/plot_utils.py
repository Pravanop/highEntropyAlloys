import matplotlib.pyplot as plt
import numpy as np
from phase_diagram_maths.math_utils import two_point_line_equation

def hmix_plotter(x_dft , y_dft , x , y, color, label, system) :
	"""

	:param x_dft:
	:param y_dft:
	:param x:
	:param y:
	"""
	
	plt.plot(x , y , c = color , alpha = 0.7 , linewidth = 4, label = "_nolegend_", zorder = 0)
	plt.scatter(x_dft , y_dft , s = 100 , c = color , edgecolors = "black", label = label)
	plt.xlabel("$x_{"+str(system[1])+"}$")
	plt.ylabel("$H_{mix}$ (meV/atom)")
	plt.title(f"{system[0]}-{system[1]}")
	plt.xlim([0,1])
	plt.ylim([0, 120])
	
	
def phase_diagram_plotter(tangent_points , system, T_range, single_energy, **kwargs) :
	"""

	:param tangent_points:
	:param system:
	"""
	xA , yA = tangent_points[: , 0 , 0] , tangent_points[: , 0 , 1]
	xB , yB = tangent_points[: , 1 , 0] , tangent_points[: , 1 , 1]
	point = None
	if 'point' in kwargs:
		point = kwargs['point']
	x_phase = np.concatenate((xA , np.flip(xB)))
	y = np.concatenate((yA , np.flip(yB)))
	max_temp = tangent_points[-1][-1][-1]
	Tm_two_point = [[0,single_energy[system[0]]['Tm']], [1, single_energy[system[1]]['Tm']]]
	x = np.linspace(0, 1, 500)
	tm_line = two_point_line_equation(two_points = Tm_two_point, x= x)
	if point:
		plt.vlines(x = point[0], ymin = T_range[0], ymax = max_temp, linestyle = "--" , colors = "#882255")
	plt.plot(x, tm_line,c = "#882255" )
	# plt.axhline(y = max_temp, linewidth = 2 , linestyle = "--" , c = "#882255")
	
	plt.scatter(x_phase , y , marker = 'o' , c = "#882255")
	plt.ylim([T_range[0] , T_range[1]])
	plt.xlim([0,1])
	plt.ylabel('T (K)')
	plt.xlabel(f'x_{system[1]}')

def gibbs_vs_T_plotter(x, gibbs_T , step, T_range) :
	"""

	:param gibbs_T:
	:param step:
	:param T_range:
	"""
	for i in range(0 , len(gibbs_T) + 1 , int(step / 75)) :
		plt.plot(x , gibbs_T[i] , linewidth = 2 , c = "#882255" , alpha = (i + 20) / 100)
	
	plt.legend(range(T_range[0] , T_range[1] , step))
	plt.xlim([0, 1])
	plt.ylim([T_range[0], T_range[1]])
	plt.ylabel("Gibbs Free Energy (meV/atom)")
	plt.xlabel("mol fraction")
	plt.axhline(y = 0)
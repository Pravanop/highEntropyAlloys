import numpy as np
from scipy.optimize import curve_fit , fsolve
from sklearn.metrics import r2_score
from phase_diagram_maths.formula import gibbs_energy , configEntropy
from itertools import combinations

def fit_parameter(y_pred , y_data) :
	"""

	:param y_pred:
	:param y_data:
	:return:
	"""
	return np.round(r2_score(y_pred = y_pred , y_true = y_data) , 2)

def solution_fit(
		func ,
		xdata ,
		ydata
		) :
	"""

	:param func:
	:param xdata:
	:param ydata:
	:return:
	"""
	popt , pcov = curve_fit(
			func , xdata = xdata , ydata = ydata
			)
	
	assert np.log10(np.linalg.cond(pcov)) < 10
	
	return fit_parameter(y_pred = func(xdata , *popt) , y_data = ydata) , popt

def get_more_points(func , points , popt) :
	"""

	:param func:
	:param points:
	:param popt:
	:return:
	"""
	xi = np.linspace(0 , 1 , points)
	return xi , func(xi , *popt)

def curve(x , y) :
	"""

	:param x:
	:param y:
	:return:
	"""
	return np.vstack((x , y)).T

def eV_to_meV(x) :
	"""

	:param x:
	:return:
	"""
	return np.round(x * 1000 , 1)

def two_point_line_equation(two_points , x) :
	"""

	:param two_points:
	:param x:
	:return:
	"""
	return ((two_points[1][1] - two_points[0][1]) / (two_points[1][0] - two_points[0][0])) * (x - two_points[0][0]) + \
		two_points[0][1]

def find_common_tangent_same(gibbs , x) :
	"""

	:param gibbs:
	:param x:
	:return:
	"""
	answer = []
	
	gibbs_curve = curve(x , gibbs)
	combs = np.array(list(combinations(gibbs_curve , 2)))
	for idx , two_point in enumerate(combs) :
		
		line = two_point_line_equation(two_points = two_point , x = x)
		diff = gibbs - line
		
		if np.any(diff < 0) :
			continue
		
		else :
			if abs(two_point[0][0] - two_point[1][0]) <= 0.1 :
				continue
			else :
				answer.append(two_point)
	
	return answer

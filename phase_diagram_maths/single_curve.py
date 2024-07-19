from phase_diagram_maths.formula import configEntropy
from phase_diagram_maths.data_utils import eV_to_meV
from phase_diagram_maths.commonTangent import CommonTangent
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

def phase_diagram(x , mix_enthalpy , T_range , single_energy, system, **kwargs) :
	"""

	:param system:
	:param single_energy:
	:param x:
	:param mix_enthalpy:
	:param T_range:
	:return:
	"""
	entropy = np.array([configEntropy([i , 1 - i]) for i in x])
	
	tangent_points = []
	gibbs_T = []
	system.reverse()
	x_out, y_out = [], []
	x = x*100
	for idx , T in enumerate(range(T_range[0] , T_range[1] , 200)) :
		
		if isinstance(mix_enthalpy, dict):
			gibbs = {}
			for key, value in mix_enthalpy.items():
				gibbs[key] = eV_to_meV(value - T*entropy)
				# gibbs[key] = gibbs[key] + 1000*single_energy[system[0]]["dft-energy"]*x + 1000*single_energy[system[1]][
				# 	"dft-energy"]*(1-x)
				
			gibbs_T.append(gibbs)
		else:
			gibbs = mix_enthalpy - T * entropy  # gibbs energy equation
			gibbs = eV_to_meV(gibbs)
			gibbs_T.append(gibbs)
			
		point = None
		if 'point' in kwargs :
			point = (kwargs['point'][0] , kwargs['point'][1] * 1000)  # eVtomeV
			ct_T = CommonTangent(
					x = x ,
					gibbs = gibbs ,
					point = point
					).get_common_tangent
		else :
			ct_T = CommonTangent(
					x = x ,
					gibbs = gibbs ,
					).get_common_tangent
		tangents = []
		for keys , value in ct_T.items() :
			if value is not None :
				for i in value :
					tangents.append(i)
		if not tangents :
			continue
		
		mega_energy = []
		for i in tangents:
			mega_energy.append(i.get_line)
		for _, value in gibbs.items():
			mega_energy.append(value)
		
		mega_energy = np.array(mega_energy)
		print(mega_energy)
		min_energy_curve = mega_energy.argmin(axis = 0)
		bounds = np.where(min_energy_curve[:-1] != min_energy_curve[1:])[0] + 1
		for l in bounds:
			x_out.append(x[l])
			y_out.append(T)
		exclude = list(set(range(0, len(min_energy_curve))) - set(bounds))
		x_points = x.copy()
		x_points[exclude] = np.nan
		y_points = [T]*len(min_energy_curve)
		tangent_points.append(list(zip(x_points, y_points)))
	return (x_out, y_out), gibbs_T

def lineLineIntersection(two_point1 , two_point2) :
	"""

	:param two_point1:
	:param two_point2:
	:return:
	"""
	# Line AB represented as a1x + b1y = c1
	a1 = two_point1[1][1] - two_point1[0][1]
	b1 = two_point1[1][0] - two_point1[0][0]
	c1 = a1 * (two_point1[0][0]) + b1 * (two_point1[0][1])
	
	# Line CD represented as a2x + b2y = c2
	a2 = two_point2[1][1] - two_point2[0][1]
	b2 = two_point2[1][0] - two_point2[0][0]
	c2 = a2 * (two_point2[0][0]) + b2 * (two_point2[0][1])
	
	determinant = a1 * b2 - a2 * b1
	
	if determinant == 0 :
		return None
	else :
		x = (b2 * c1 - b1 * c2) / determinant
		y = (a1 * c2 - a2 * c1) / determinant
		return x , y

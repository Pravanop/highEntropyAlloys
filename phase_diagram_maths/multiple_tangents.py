import numpy as np
from phase_diagram_maths.commonTangent import curve , two_point_line_equation
from phase_diagram_maths.formula import (
	binary_regular_model_enthalpy , binary_subregular_model_enthalpy ,
	intermetallic_energy , kopp_neumann_law ,
	)
from phase_diagram_maths.plot_utils import hmix_plotter , phase_diagram_plotter , gibbs_vs_T_plotter
from phase_diagram_maths.data_utils import get_enthalpy_fit , process_inputs , eV_to_meV
from phase_diagram_maths.single_curve import phase_diagram_oneCurve
import matplotlib.pyplot as plt
from phase_diagram_maths.formula import get_mixing_enthalpy_system , configEntropy
from phase_diagram_maths.data_utils import solution_fit , eV_to_meV , get_enthalpy_fit
from phase_diagram_maths.commonTangent import find_common_tangent_same, find_common_tangent_pointCurve, two_point_line_equation
from vasp_input_main.file_utils import load_json_to_dict

alloy = "Cr-Ta"
lattice = "BCC"
mol_fraction , dft_energies , system = process_inputs(alloy = alloy , lattice = lattice)
single_energy_path = "/Users/pravanomprakash/Documents/Projects/highEntropyAlloys/data/single_energy.json"
single_energy = load_json_to_dict(single_energy_path)
alloy_enthalpies , x , mix_enthalpy = get_enthalpy_fit(
		mol_fraction = mol_fraction , dft_energies = dft_energies , system = system , no_atoms = 24 ,
		lattice = lattice , model = binary_regular_model_enthalpy, single_energy= single_energy
		)

def phase_diagram_pointCurve(x , mix_enthalpy , T_range , point) :
	"""

	:param x:
	:param mix_enthalpy:
	:param T_range:
	:return:
	"""
	entropy = np.array([configEntropy([i , 1 - i]) for i in x])
	
	tangent_points = []
	gibbs_T = []
	im_T = []
	for idx , T in enumerate(range(T_range[0] , T_range[1] , 75)) :
		gibbs = mix_enthalpy - T * entropy  # gibbs energy equation
		# im = eV_to_meV(intermetallic_energy(a = Cp_im[0], b = Cp_im[1], c = Cp_im[2], T = T, H0 = point[1]))
		im = eV_to_meV(point[1])
		im_T.append(im)
		gibbs = eV_to_meV(gibbs)
		gibbs_T.append(gibbs)
		ct_pointT = find_common_tangent_pointCurve(gibbs , x , [point[0], im])
		ct_T = find_common_tangent_same(gibbs, x)
		
		break
		# plt.plot(x,gibbs)
		# plt.plot(x,linePoint)
		# plt.plot(x,lineT)
		# plt.show()
	return np.array(tangent_points) , np.array(gibbs_T), np.array(im_T)


point = [0.33 , -0.09]
mass = [single_energy[i]['mass'] for i in system]
cp = np.array([[single_energy[i]['a'], single_energy[i]['b'], single_energy[i]['c']] for i in system])
intermetallic_cp = kopp_neumann_law(massA = mass[0], massB = mass[1], cpA = cp[0], cpB = cp[1], mol_fraction = point[0])
# print(intermetallic_cp, cp)
T_range = [500 , 3000]
common_tangent , gibbs_T, im_T = phase_diagram_pointCurve(x , mix_enthalpy , T_range , point, intermetallic_cp)
# print(common_tangent)
phase_diagram_plotter(common_tangent, system, T_range)
plt.show()
# for i in range(len(im_T)):
# 	plt.scatter(point[0] , im_T[i])
# gibbs_vs_T_plotter(x , gibbs_T , 300 , T_range)
# plt.show()

from phase_diagram_maths.formula import get_mixing_enthalpy_system , configEntropy
from phase_diagram_maths.data_utils import solution_fit , eV_to_meV , find_common_tangent_same
from formula import binary_regular_model_enthalpy, binary_subregular_model_enthalpy
from plot_utils import hmix_plotter, phase_diagram_plotter, gibbs_vs_T_plotter
import numpy as np
import matplotlib.pyplot as plt
from vasp_input_main.file_utils import load_json_to_dict
import pandas as pd

def get_enthalpy_fit(mol_fraction , single_energy_path , dft_energies , system , no_atoms , lattice , model) :
	"""

	:param mol_fraction:
	:param single_energy_path:
	:param dft_energies:
	:param system:
	:param no_atoms:
	:param lattice:
	:param model:
	:return:
	"""
	single_energy = load_json_to_dict(single_energy_path)
	
	end_member_energy = [single_energy[i] for i in system]
	alloy_enthalpies = get_mixing_enthalpy_system(
			dft_alloy_energies = dft_energies ,
			mol_fraction = mol_fraction ,
			system = system ,
			end_member_energy = end_member_energy ,
			no_atoms = no_atoms ,
			lattice = lattice
			)
	
	alloy_enthalpies = np.pad(alloy_enthalpies , (1 , 1) , 'constant')
	r2 , popt = solution_fit(xdata = mol_fraction , ydata = alloy_enthalpies , func = model)
	print("R2:" , r2)
	
	x = np.linspace(0 , 1 , 101)
	mix_enthalpy = model(x , *popt)
	
	return alloy_enthalpies , x , mix_enthalpy

def phase_diagram_oneCurve(x , mix_enthalpy , T_range) :
	"""

	:param x:
	:param mix_enthalpy:
	:param T_range:
	:return:
	"""
	entropy = np.array([configEntropy([i , 1 - i]) for i in x])
	
	tangent_points = []
	gibbs_T = []
	for idx , T in enumerate(range(T_range[0] , T_range[1] , 75)) :
		
		gibbs = mix_enthalpy - T * entropy #gibbs energy equation
		gibbs = eV_to_meV(gibbs)
		gibbs_T.append(gibbs)
		ct_T = find_common_tangent_same(gibbs , x)
		if ct_T :
			x1 , x2 = ct_T[0][0][0] , ct_T[0][1][0]
			T1 , T2 = T , T
			tangent_points.append([[x1 , T1] , [x2 , T2]])
	
	return np.array(tangent_points) , np.array(gibbs_T)


alloy = 'Cr-W'
lattice = "BCC"
data_path = f"/Users/pravanomprakash/Documents/Projects/highEntropyAlloys/data/binary_{lattice}.csv"
df_BCC = pd.read_csv(
	data_path , header = 0 ,
	index_col = "Unnamed: 0"
	)
mol_fraction = np.append(np.array(df_BCC.columns).astype(float) , 1)
mol_fraction = np.append(0 , mol_fraction)
dft_energies = np.array(df_BCC.loc[alloy].to_list())

system = alloy.split('-')
system.reverse()
single_energy_path = "/Users/pravanomprakash/Documents/Projects/highEntropyAlloys/data/single_energy.json"

alloy_enthalpies , x , mix_enthalpy = get_enthalpy_fit(
		mol_fraction = mol_fraction , single_energy_path =
		single_energy_path , dft_energies = dft_energies , system = system , no_atoms = 24 , lattice = lattice ,
		model = binary_regular_model_enthalpy
		)

hmix_plotter(mol_fraction , eV_to_meV(alloy_enthalpies) , x , eV_to_meV(mix_enthalpy))
plt.show()
T_range = (500 , 3000)
tangent_points , gibbs_T = phase_diagram_oneCurve(x , mix_enthalpy , T_range = T_range)
phase_diagram_plotter(tangent_points , system)
plt.show()
gibbs_vs_T_plotter(gibbs_T , 300)
plt.show()

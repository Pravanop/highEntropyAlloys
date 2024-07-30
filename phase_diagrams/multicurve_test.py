import numpy as np
import itertools
from phase_diagram_maths.formula import binary_regular_model_enthalpy , binary_subregular_model_enthalpy
from phase_diagram_maths.plot_utils import hmix_plotter , phase_diagram_plotter , gibbs_vs_T_plotter
from phase_diagram_maths.data_utils import get_enthalpy_fit , process_inputs, eV_to_meV, load_json_to_dict, flatten
import matplotlib.pyplot as plt
from phase_diagram_maths.single_curve import phase_diagram

single_energy = load_json_to_dict()
alloy = "Cr-W"
lattice = "BCC"
alloy2 = "Cr-Ta"
lattice2 = "BCC"
mol_fraction , dft_energies , system = process_inputs(alloy = alloy , lattice = lattice)
mol_fraction2 , dft_energies2 , system2 = process_inputs(alloy = alloy2 , lattice = lattice2)

alloy_enthalpies1 , x , mix_enthalpy1 = get_enthalpy_fit(
		mol_fraction = mol_fraction ,
		dft_energies = dft_energies ,
		system = system ,
		no_atoms = 24 ,
		lattice = lattice ,
		model = binary_subregular_model_enthalpy,
		single_energy = single_energy
		)
alloy_enthalpies2 , _ , mix_enthalpy2 = get_enthalpy_fit(
		mol_fraction = mol_fraction2 ,
		dft_energies = dft_energies2 ,
		system = system2 ,
		no_atoms = 24 ,
		lattice = lattice2 ,
		model = binary_subregular_model_enthalpy,
		single_energy = single_energy
		)

mix_enthalpy = {
		"BCC":mix_enthalpy1,
		"HCP":mix_enthalpy2[::-1]
		}
point = (0.5, -0.05)
T_range = (500 , max(single_energy[system[0]]['Tm'], single_energy[system[1]]['Tm'])+50)
tangent_points = np.array(phase_diagram(x , mix_enthalpy , T_range = T_range, point = point))
new_pad = []
for idx , i in enumerate(tangent_points) :

	new_pad.append(i.T)

new_pad = np.array(new_pad)
# print(new_pad)
x , y = [] , []
for idx , i in enumerate(new_pad) :
	# print(new_pad[idx][0])
	x.append(new_pad[idx][0])
	y.append(new_pad[idx][1])

x = np.array(x).T
y = np.array(y).T
# print(x)
# print(y)
# print(x)
for i in range(len(x)) :
	plt.scatter(x[i] , y[i], c= 'black')
# plt.imshow(y)
plt.xlim([0,1])
plt.ylim([T_range[0],T_range[1]])
plt.xlabel(f'$x_{system[1]}$')
plt.ylabel("T (K)")
plt.show()
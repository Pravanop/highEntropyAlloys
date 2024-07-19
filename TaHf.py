import itertools

from phase_diagram_maths.formula import binary_regular_model_enthalpy , binary_subregular_model_enthalpy
from phase_diagram_maths.plot_utils import hmix_plotter , phase_diagram_plotter , gibbs_vs_T_plotter
from phase_diagram_maths.data_utils import get_enthalpy_fit , process_inputs , eV_to_meV , load_json_to_dict
from phase_diagram_maths.single_curve import phase_diagram
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

font = {
		'family' : 'Helvetica' ,
		'size'   : 16
		}
matplotlib.rc('font' , **font)
single_energy = load_json_to_dict()
alloy = "Ta-Hf"
lattice = ["BCC" , "HCP"]
color = {"BCC": "blue",
         "HCP":"red"}
mol_fraction = {}
dft_energies = {}
for i in lattice:
	mol_fraction[i] , dft_energies[i] , system = process_inputs(alloy = alloy , lattice = i)

mix_enthalpy = {}
alloy_enthalpy = {}
for i in lattice :
	alloy_enthalpy[i] , x , mix_enthalpy[i] = get_enthalpy_fit(
			mol_fraction = mol_fraction[i] ,
			dft_energies = dft_energies[i] ,
			system = system ,
			no_atoms = 24 ,
			lattice = lattice ,
			model = binary_subregular_model_enthalpy ,
			single_energy = single_energy
			)
print(alloy_enthalpy)
for key , value in mix_enthalpy.items() :
	hmix_plotter(
		mol_fraction[key] , eV_to_meV(alloy_enthalpy[key]) , x , eV_to_meV(value) , color = color[key] , label = key ,
		system = system
		)
plt.legend(lattice)
# plt.show()
T_range = (500 , max(single_energy[system[0]]['Tm'] , single_energy[system[1]]['Tm']) + 50)
tangent_points, gibbs = phase_diagram(x , mix_enthalpy , T_range = T_range, single_energy = single_energy, system = system)
gibbs_plot = gibbs[0::1]
fig, axs = plt.subplots(1, 1, sharex = True, sharey = True, figsize = (10, 8))
c = plt.cm.Blues(np.linspace(0.5,1, len(gibbs_plot)))
for idx, i in enumerate(gibbs_plot):
	temp = i["BCC"]
	axs.plot(x, temp, color = c[idx], linewidth =  4)
	axs.set_xlim([0, 1])
sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=T_range[0], vmax=T_range[1]))
cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
cbar = fig.colorbar(sm, cax=cbar_ax, aspect = 0.3)
cbar.ax.set_ylabel("Temperature")
axs.set_xlabel("$X_{Cr}$")
axs.set_ylabel("G (meV/atom)")
rang = abs(T_range[0] - T_range[1])
# plt.show()
tangent_points = np.array(tangent_points)
new_pad = []
for idx , i in enumerate(tangent_points) :
	
	new_pad.append(i.T)

new_pad = np.array(new_pad)

x , y = [] , []
for idx , i in enumerate(new_pad) :
	x.append(new_pad[idx][0])
	y.append(new_pad[idx][1])

x = np.array(x).T
y = np.array(y).T
print(x)
for i in range(len(x)) :
	plt.scatter(x[i] , y[i] , c = 'black')
# plt.imshow(y)
# plt.xlim([0 , 1])
# plt.ylim([T_range[0] , T_range[1]])
# plt.xlabel(f'x_{system[1]}')
# plt.ylabel("T (K)")
# plt.show()

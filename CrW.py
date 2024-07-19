import numpy as np

from phase_diagram_maths.formula import binary_regular_model_enthalpy , binary_subregular_model_enthalpy
from phase_diagram_maths.plot_utils import hmix_plotter , phase_diagram_plotter , gibbs_vs_T_plotter
from phase_diagram_maths.data_utils import get_enthalpy_fit , process_inputs , eV_to_meV , load_json_to_dict, solution_fit
from phase_diagram_maths.single_curve import phase_diagram
import matplotlib.pyplot as plt
import matplotlib
# from csaps import csaps

font = {
		'family' : 'Helvetica' ,
		'size'   : 20
		}
matplotlib.rc('font' , **font)

alloy = "Cr-W"
lattice = "BCC"
mol_fraction , dft_energies , system = process_inputs(alloy = alloy , lattice = lattice)
single_energy_path = "/Users/pravanomprakash/Documents/Projects/highEntropyAlloys/data/single_energy.json"
single_energy = load_json_to_dict(single_energy_path)

alloy_enthalpies , x , mix_enthalpy = get_enthalpy_fit(
		mol_fraction = mol_fraction , dft_energies = dft_energies , system = system , no_atoms = 24 ,
		lattice = lattice ,
		model = binary_regular_model_enthalpy , single_energy = single_energy
		)

# Part 1
mix_enthalpy = {
		"BCC" : mix_enthalpy ,
		}
# for key , value in mix_enthalpy.items() :
# 	hmix_plotter(
# 		mol_fraction , eV_to_meV(alloy_enthalpies) , x , eV_to_meV(value) , color = '#CC6677' , label = key ,
# 		system = system
# 		)
# plt.show()

T_range = (500 , 3000 + 50)
tangent_points , gibbs = phase_diagram(x , mix_enthalpy , T_range = T_range, single_energy = single_energy, system = system)
gibbs_plot = gibbs[0::1]
# fig, axs = plt.subplots(1, 1, sharex = True, sharey = True, figsize = (10, 8))
# c = plt.cm.Blues(np.linspace(0.5,1, len(gibbs_plot)))
# for idx, i in enumerate(gibbs_plot):
# 	temp = i["BCC"]
# 	axs.plot(x, temp, color = c[idx], linewidth =  4)
# 	axs.set_xlim([0, 1])
# sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=T_range[0], vmax=T_range[1]))
# cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
# cbar = fig.colorbar(sm, cax=cbar_ax, aspect = 0.3)
# cbar.ax.set_ylabel("Temperature")
# axs.set_xlabel("$X_{Cr}$")
# axs.set_ylabel("G (meV/atom)")
# rang = abs(T_range[0] - T_range[1])
# plt.show()
print(tangent_points)
xy = list(zip(tangent_points[0], tangent_points[1]))
# xy.sort()
xy = np.array(xy)
sort = xy[np.lexsort((xy[:,1],xy[:,0]))]
print(sort)
x, y = sort[:,0], sort[:, 1]
plt.plot(x,y, linewidth = 4)
plt.xlim([0,100])
plt.ylim([1000, T_range[1]])
plt.ylabel("T(K)")
plt.xlabel("$x_{W}$")
plt.show()
# tangent_points = np.array(tangent_points)
# new_pad = []
# for idx , i in enumerate(tangent_points) :
#
# 	new_pad.append(i.T)
#
# new_pad = np.array(new_pad)
#
# x , y = [] , []
# for idx , i in enumerate(new_pad) :
# 	# print(new_pad[idx][0])
# 	x.append(new_pad[idx][0])
# 	y.append(new_pad[idx][1])
#
# x = np.array(x).T
# y = np.array(y).T
# print(x)
# for i in range(len(x)) :
# 	plt.scatter(x[i] , y[i] , c = 'black')
# # plt.imshow(y)
# plt.xlim([0 , 1])
# plt.ylim([T_range[0] , T_range[1]])
# plt.xlabel(f'x_{system[1]}')
# plt.ylabel("T (K)")
# plt.show()

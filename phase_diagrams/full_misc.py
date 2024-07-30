import numpy as np

from phase_diagram_maths.formula import binary_regular_model_enthalpy , binary_subregular_model_enthalpy
from phase_diagram_maths.plot_utils import hmix_plotter , phase_diagram_plotter , gibbs_vs_T_plotter
from phase_diagram_maths.data_utils import get_enthalpy_fit , process_inputs , eV_to_meV , load_json_to_dict
from phase_diagram_maths.single_curve import phase_diagram
import matplotlib.pyplot as plt
import matplotlib

font = {
		'family' : 'Helvetica' ,
		'size'   : 20
		}
matplotlib.rc('font' , **font)

# alloys = ["Cr-V" , "Ta-W" , "V-W"]
alloys = ["Cr-W" , "Cr-Ta" , "V-W"]
lattice = "BCC"

mix_enthalpy_dict = {}
for alloy in alloys :
	print(alloy)
	mol_fraction , dft_energies , system = process_inputs(alloy = alloy , lattice = lattice)
	single_energy_path = "/data/single_energy.json"
	single_energy = load_json_to_dict(single_energy_path)
	
	alloy_enthalpies , x , mix_enthalpy = get_enthalpy_fit(
			mol_fraction = mol_fraction , dft_energies = dft_energies , system = system , no_atoms = 24 ,
			lattice = lattice ,
			model = binary_regular_model_enthalpy , single_energy = single_energy
			)
	mix_enthalpy_dict[alloy] = {
			"alloy_enthalpies" : eV_to_meV(alloy_enthalpies) ,
			"mix_enthalpy"     : eV_to_meV(mix_enthalpy) ,
			"system"           : system
			}
c = plt.cm.coolwarm(np.linspace(0.0 , 1 , 3))
count = 0
fig, axs = plt.subplots(nrows = 3, ncols = 1, sharex = True, sharey = True, figsize = (10 , 8))
for key , value in mix_enthalpy_dict.items() :
	print(key)
	axs[count].plot(
		x , value['mix_enthalpy'] , color = c[count] , alpha = 0.7 , linewidth = 4 , label = "_nolegend_" ,
		zorder = 0
		)
	value['system'].reverse()
	axs[count].scatter(mol_fraction , value['alloy_enthalpies'] , s = 100 , c = c[count] , edgecolors = "black")
	axs[count].text(x = 0.15, y = 20, s = '-'.join(value['system']))
	axs[count].set_xlabel("$x_{" + str(value['system'][1]) + "}$")
	axs[count].set_xlim([0 , 1])
	axs[count].set_ylim([-100 , 120])
	axs[count].axhline(y = 0, linestyle = "-", linewidth = 1, color = "black")
	count += 1

plt.subplots_adjust(left = 0.13, bottom = 0.10, right = 0.47, top = 0.84, wspace = 0.22, hspace = 0.23)
fig.supylabel("$H_{mix}$ (meV/atom)")
# plt.title(f"{value['system'][0]}-{value['system'][1]}")
plt.savefig("/Users/pravanomprakash/Documents/Projects/highEntropyAlloys/Monte Carlo/plots/example_energy.png", dpi = 300)

plt.show()

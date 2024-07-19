import json

import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from itertools import combinations
from phase_diagram_maths.formula import get_mixing_enthalpy_system
import pandas as pd

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


def eV_to_meV(x) :
	"""

	:param x:
	:return:
	"""
	return np.round(x * 1000 , 2)


def get_enthalpy_fit(
		mol_fraction , dft_energies , system ,single_energy, no_atoms , lattice , model
		) :
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
	
	end_member_energy = [single_energy[i]["dft-energy"] for i in system]
	alloy_enthalpies = get_mixing_enthalpy_system(
			dft_alloy_energies = dft_energies ,
			mol_fraction = mol_fraction ,
			system = system ,
			end_member_energy = end_member_energy ,
			no_atoms = no_atoms ,
			lattice = lattice
			)
	
	alloy_enthalpies = np.pad(alloy_enthalpies , (1 , 1) , 'constant')
	ydata = alloy_enthalpies.copy()
	isnan = np.isnan(ydata)
	ydata = ydata[~isnan]
	xdata = mol_fraction.copy()
	xdata = xdata[~isnan]
	r2 , popt = solution_fit(xdata = xdata , ydata = ydata , func = model)
	print("R2:" , r2)
	
	x = np.linspace(0 , 1 , 71)
	# print("Omega1: ",popt[0], " Omega2: ", popt[1])
	mix_enthalpy = model(x , *popt)
	
	return alloy_enthalpies , x , mix_enthalpy

def process_inputs(
		alloy: str , lattice: str , data_path: str =
		"/Users/pravanomprakash/Documents/Projects/highEntropyAlloys/data"
		) :
	data_path = f"{data_path}/binary_{lattice}.csv"
	df_BCC = pd.read_csv(
			data_path , header = 0 ,
			index_col = "Unnamed: 0"
			)
	mol_fraction = np.append(np.array(df_BCC.columns).astype(float) , 1)
	mol_fraction = np.append(0 , mol_fraction)
	dft_energies = np.array(df_BCC.loc[alloy].to_list())
	
	system = alloy.split('-')
	system.reverse()
	
	return mol_fraction , dft_energies , system

def load_json_to_dict(json_file_path = "/Users/pravanomprakash/Documents/Projects/highEntropyAlloys/data/single_energy.json"):
    try:
        with open(json_file_path, 'r') as file:
            data = json.load(file)
        return data
    except json.JSONDecodeError:
        print("Error: The file is not a valid JSON.")
    except FileNotFoundError:
        print(f"Error: The file {json_file_path} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def flatten(xss) :
	"""

	:param xss:
	:return:
	"""
	return [x for xs in xss for x in xs]
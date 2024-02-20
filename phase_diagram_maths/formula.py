import numpy as np

def mixing_enthalpy(
		alloy_energy: float ,
		mol_fraction: np.array ,
		ele_list: list[str] ,
		end_member_energy: np.array
		) -> float :
	"""
	Calculates mixing enthalpy of alloys with DFT obtained energies.
					delH_mix = H_alloy - sum_i^n(x_i*H_i)
	where x is mol fractions, and all the RHS energies are DFT obtained values.
	:param alloy_energy: The DFT obtained alloy energy in ev/atom
	:param mol_fraction: A list of float values whose sum should equal one
	:param ele_list: A list of strings with element names, just for reference
	:param end_member_energy: A list of single element energies in ev/atom, must be of same order as mol_fraction
	:return: mixing enthalpy value in ev/atom
	"""
	return np.round(alloy_energy - np.sum(mol_fraction * end_member_energy) , 3)

def batch_mixing_enthalpy_bycomp(
		alloy_energies: np.array ,
		mol_fraction: np.array(np.array) ,
		ele_list: list[str] ,
		end_member_energy: np.array
		) -> np.array :
	"""
	For a certain system, mixing enthalpy is calculated across a range of composition.
	:param alloy_energies: The set of DFT obtained alloy energies in ev/atom
	:param mol_fraction: The set of mol_fraction in the same order as DFT alloys
	:param ele_list: A list of strings with element names, just for reference
	:param end_member_energy: A list of single element energies in ev/atom, must be of same order as mol_fraction
	:return: an array of mixing enthalpy values in ev/atom
	"""
	alloy_enthalpies = []
	for idx, alloy_energy in enumerate(alloy_energies):
		alloy_enthalpies.append(mixing_enthalpy(alloy_energy, np.array([mol_fraction[idx+1], 1- mol_fraction[
			idx+1]]), ele_list = ele_list, end_member_energy = end_member_energy))
	return np.array(alloy_enthalpies)

def gibbs_energy(
		enthalpy: float ,
		entropy: float ,
		temperature: float
		) -> float :
	"""
	Function to implement G = H - TS
	:param enthalpy: the mixing enthalpy
	:param entropy: the configurational entropy
	:param temperature: the desired temperature
	:return: gibbs free energy
	"""
	return np.round(enthalpy - temperature * entropy , 7)

def configEntropy(mol_ratio: np.array) -> float :
	"""
	A simple function to calculate boltzmann configurational entropy.

					delta_S = -k_b*sum_i=1_n(x_i*ln(x_i))

	:param mol_ratio: list of floats of mole fraction
	:return: boltzmann configurational entropy
	"""
	
	k_b = 8.617333262e-05
	if 1.0 in mol_ratio :
		return 0
	return np.round(-k_b * np.sum(mol_ratio * np.log(mol_ratio)) , 7)

def binary_regular_model_enthalpy(
		mol_fraction: np.array ,
		omega: float
		) -> np.array :
	"""
	Basic function to calculate a regular model. Taken from John Cavin's thesis Eq 3.2.
	:param omega: array of pairwise-interaction parameters for each mole fraction
	:param mol_fraction: array of mole fractions
	:return: mixing enthalpy value for mole fraction
	"""
	x_i = mol_fraction
	return x_i * (1 - x_i) * omega

def binary_subregular_model_enthalpy(
		mol_fraction: np.array ,
		omega1: float ,
		omega2: float
		) -> np.array :
	"""
	Basic function to calculate a subregular model with cubic fit. Taken from John Cavin's thesis Eq 3.3.
	:param omega: array of pairwise-interaction parameters for each mole fraction
	:param mol_fraction: array of mole fractions
	:return: mixing enthalpy value for mole fraction
	"""
	x_i = mol_fraction
	return x_i * (1 - x_i) * (omega1 * x_i + omega2 * (1 - x_i))

def binary_subregular4_model_enthalpy(
		mol_fraction: np.array ,
		omega1: float ,
		omega2: float ,
		omega3: float
		) -> np.array :
	"""
	Basic function to calculate a subregular model with cubic fit. Taken from John Cavin's thesis Eq 3.4.
	:param omega: array of pairwise-interaction parameters for each mole fraction
	:param mol_fraction: array of mole fractions
	:return: mixing enthalpy value for mole fraction
	"""
	x_i = mol_fraction
	return x_i * (1 - x_i) * (omega1 * (1 - x_i) + omega2 * (np.power((1 - x_i) , 2)) + omega3)

def dft_energy_per_atom(
		alloy_energies: np.array ,
		no_atoms: int
		) -> np.array :
	
	"""
	
	:param alloy_energies:
	:param no_atoms:
	:return:
	"""
	
	return alloy_energies / no_atoms

def get_mixing_enthalpy_system(
		dft_alloy_energies: np.array ,
		mol_fraction: np.array ,
		system: list[str] ,
		end_member_energy: np.array ,
		no_atoms: int ,
		lattice: str
		) -> np.array :
	
	"""
	
	:param lattice:
	:param no_atoms:
	:param dft_alloy_energies:
	:param mol_fraction:
	:param system:
	:param end_member_energy:
	:return:
	"""
	dft_energy_atom = dft_energy_per_atom(dft_alloy_energies , no_atoms)
	return batch_mixing_enthalpy_bycomp(
		alloy_energies = dft_energy_atom ,
		mol_fraction = mol_fraction ,
		ele_list = system ,
		end_member_energy = end_member_energy
		)

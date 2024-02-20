import numpy as np

def deriv_binaryconfigEntropy(mol_ratio: np.array) -> float :
	"""

	:param mol_ratio:
	:return:
	"""
	k_b = 8.617333262e-05
	
	if 1.0 in mol_ratio :
		return 0
	x_i = mol_ratio[0]
	return np.round(-k_b * (np.log(x_i) - (np.log(1 - x_i))) , 7)

def deriv_binary_subregular4_model_enthalpy(
		mol_fraction: np.array ,
		omega1: float ,
		omega2: float ,
		omega3: float
		) -> np.array :
	x_i = mol_fraction
	return x_i * (1 - x_i) * (-omega1 - 2 * omega2 * (1 - x_i)) + (
			omega1 * (1 - x_i) + omega2 * (np.power((1 - x_i) , 2)) + omega3) * (1 - 2 * x_i)

def deriv_binary_subregular_model_enthalpy(
		mol_fraction: np.array ,
		omega1: float ,
		omega2: float
		) -> np.array :
	"""
	Derivative for a binary subregular model with cubic fit.
	:param mol_fraction:
	:param omega1:
	:param omega2:
	"""
	x_i = mol_fraction
	return x_i * (1 - x_i) * (omega1 - omega2) + (1 - 2 * x_i) * (omega1 * x_i + omega2 * (1 - x_i))

def deriv_binary_regular_model_enthalpy(
		mol_fraction: np.array ,
		omega: float
		) -> np.array :
	"""
	Derivative for a binary regular model.
	:param mol_fraction:
	:param omega:
	:return:
	"""
	x_i = mol_fraction
	return (1 - 2 * x_i) * omega
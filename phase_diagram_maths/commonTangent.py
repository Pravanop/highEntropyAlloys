import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from itertools import combinations

def gradient(y, x, step):
	dx = np.max(x) - np.min(x)
	dx = dx/step
	return np.round(np.gradient(y, dx), 2)

def solve_common_tangent_same(curve,x, step):
	grad = gradient(curve, x, step)
	grad_copy = grad
	curve = np.round(curve, 4)
	for idx, i in enumerate(grad):
		tangent = line_equation(i, (x[idx],curve[idx]), x)
		inter = set(curve).intersection(set(tangent))
		if len(list(inter)) == 1:
			print(x[curve == float(list(inter)[0])])

def line_equation(slope, points, x):
	return np.round(slope*(x - points[0]) + points[1],4)


def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')


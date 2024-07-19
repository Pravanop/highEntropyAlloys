import numpy as np
import warnings
def movingaverage(interval , window_size) :
	"""

    :param interval:
    :param window_size:
    :return:
    """
	window = np.ones(int(window_size)) / float(window_size)
	return np.convolve(interval , window , 'same')

def curve(x , y) -> np.array :
	"""

	:param x:
	:param y:
	:return: a numpy array of [x,y] points
	"""
	return np.vstack((x , y)).T

def slope_line_equation(slope , points , x) :
	"""

	:param slope:
	:param points:
	:param x:
	:return:
	"""
	return np.round(slope * (x - points[0]) + points[1] , 4)

def two_point_line_equation(two_points: np.array , x) :
	"""

	:param two_points: A 2,2 matrix of form [[x1, y1], [x2, y2]]
	:param x: the range of x points
	:return:
	"""
	warnings.filterwarnings("ignore")
	return ((two_points[1][1] - two_points[0][1]) / (two_points[1][0] - two_points[0][0])) * (x - two_points[0][0]) + \
		two_points[0][1]

def gradient(y , x , step) :
	dx = np.max(x) - np.min(x)
	dx = dx / step
	return np.round(np.gradient(y , dx) , 2)

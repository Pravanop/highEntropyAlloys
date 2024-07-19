import numpy as np

class tangent :
	"""
	Tangent class that takes in two points and a range of mol fraction to create a line.
	"""
	
	def __init__(self , two_point , line , x) :
		self.two_point = two_point
		self.line = line
		self.x = x
	
	@property
	def get_points(self) -> np.array :
		"""
		Get the two points at a later time
		:return: 2D array of [[x1, y1],[x2,y2]]
		"""
		return self.two_point
	
	@property
	def get_line(self) -> np.array :
		"""
		get the line at a later time
		:return: 1D array of points for a given x
		"""
		x1 , x2 = min(self.two_point[0][0] , self.two_point[1][0]) , max(self.two_point[0][0] , self.two_point[1][0])
		line_trim = self.line
		line_trim[(self.x < x1) | (self.x > x2)] = np.inf
		return line_trim

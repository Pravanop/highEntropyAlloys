import numpy as np
from itertools import combinations , product
from phase_diagram_maths.math_utils import curve , two_point_line_equation
from phase_diagram_maths.tangent import tangent

class CommonTangent :
	
	def __init__(self , x , gibbs , **kwargs) :
		
		self.intermetallic = None
		self.intermetallic_list = None
		if "point" in kwargs :
			if isinstance(kwargs["point"] , list) :
				self.intermetallic_list = kwargs['point']
			else :
				self.intermetallic = kwargs["point"]
		
		self.gibbs = gibbs
		self.x = x
		self.tolerance = 7
	
	def find_common_tangent_intermetallic(self , gibbs , point) :
		"""

		:param gibbs:
		:param point:
		:return:
		"""
		answer = []
		gibbs_curve = curve(self.x , gibbs)
		for idx , curve_point in enumerate(gibbs_curve) :
			if curve_point[0] < 1 or curve_point[0] > 99:
				continue
			two_point = np.hstack((curve_point , point)).reshape(2 , 2)
			line = two_point_line_equation(two_points = two_point , x = self.x)
			diff = gibbs - line
			print(diff)
			if np.any(diff < 0) :
				continue
				
			else :
				answer.append(tangent(two_point , line , self.x))
		
		if not answer :
			return None
		else :
			return answer
	
	def find_common_tangent_same(self , gibbs) :
		"""

		:param gibbs:
		:param tolerance:
		:return:
		"""
		answer = []
		gibbs_curve = curve(self.x , gibbs)
		combs = np.array(list(combinations(gibbs_curve , 2)))
		for idx , two_point in enumerate(combs) :
			if two_point[0][0] < 0 or two_point[0][0] > 100:
				continue
			if two_point[1][0] < 0 or two_point[1][0] > 100:
				continue
			line = two_point_line_equation(two_points = two_point , x = self.x)
			diff = gibbs - line
			# print(diff)
			if np.any(diff < 0) :
				continue
			
			else :
				if abs(two_point[0][0] - two_point[1][0]) <= self.tolerance :
					# answer.append(tangent(two_point , line , self.x))
					continue
				else :
					answer.append(tangent(two_point , line , self.x))
		
		if not answer :
			return None
		else :
			return answer
	
	def find_common_tangent_intermetallic_list(self , gibbs) :
		"""

		:return:
		"""
		intermetallic_tangent_list = []
		for idx , i in enumerate(self.intermetallic_list) :
			intermetallic_tangent_list.append(self.find_common_tangent_intermetallic(gibbs , i))
		
		return intermetallic_tangent_list
	
	def find_common_tangent_twocurve(self , gibbs_pair) :
		"""

		:param gibbs_pair:
		:param tolerance:
		:return:
		"""
		answer = []
		combs = np.array((list(product(curve(self.x , gibbs_pair[0]) , curve(self.x , gibbs_pair[1])))))
		for idx , two_point in enumerate(combs) :
			if two_point[0][0] < 0.01 or two_point[0][0] > 0.99 :
				continue
			if two_point[1][0] < 0.01 or two_point[1][0] > 0.99 :
				continue
			line = two_point_line_equation(two_points = two_point , x = self.x)
			diff1 = gibbs_pair[0] - line
			diff2 = gibbs_pair[1] - line
			
			if np.any(diff1 < 0) or np.any(diff2 < 0) :
				continue
			
			else :
				if abs(two_point[0][0] - two_point[1][0]) <= self.tolerance :
					# answer.append(tangent(two_point , line , self.x))
					continue
				else :
					answer.append(tangent(two_point , line , self.x))
		
		if not answer :
			return None
		else :
			return answer
	
	def find_common_tangent_multicurve(self) :
		key_combs = ['-'.join(i) for i in list(combinations(list(self.gibbs.keys()) , 2))]
		key_combs = key_combs + list(self.gibbs.keys())
		tangent_dict = {}
		for idx , inst in enumerate(key_combs) :
			if '-' in inst :
				keys = inst.split('-')
				gibbs_pair = (self.gibbs[keys[0]] , self.gibbs[keys[1]])
				temp = self.find_common_tangent_twocurve(gibbs_pair = gibbs_pair)
				tangent_dict[inst] = temp
			else :
				# print(self.gibbs[inst])
				temp = self.find_common_tangent_same(self.gibbs[inst])
				tangent_dict[inst] = temp
				if self.intermetallic :
					tangent_dict["intermetallic"] = self.find_common_tangent_intermetallic(
							self.gibbs[inst] ,
							self.intermetallic
							)
				
				if self.intermetallic_list :
					tangent_dict["intermetallic"] = self.find_common_tangent_intermetallic_list()
		
		return tangent_dict
	
	@property
	def get_common_tangent(self) -> dict :
		"""

		:return:
		"""
		if not isinstance(self.gibbs , dict) :
			tangent_dict = {}
			
			if self.intermetallic :
				tangent_dict["intermetallic"] = self.find_common_tangent_intermetallic(self.intermetallic)
			
			if self.intermetallic_list :
				tangent_dict["intermetallic"] = self.find_common_tangent_intermetallic_list()
			
			tangent_dict["same"] = [self.find_common_tangent_same()]
		else :
			tangent_dict = self.find_common_tangent_multicurve()
		return tangent_dict

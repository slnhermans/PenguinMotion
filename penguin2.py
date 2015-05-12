import numpy as np
import scipy as sp
import random

class Penguin(object):

	def __init__(self, radius, position, alignment, boundary):
		"""
		Arguments:
			radius		:=	float representing penguin radius
			position	:=	numpy array representing penguin position (x, y, z)
			alignment	:=	numpy array representing alignment unit vector
			boundary	:=	boolean value, True if boundary penguin, False if bulk penguin
		"""

		self.radius = radius
		self.position = position
		self.alignment = alignment
		self.boundary = boundary

	def get_distance(self, penguin2):

		return self.position - penguin2.position

	def net_force(self, F_self, F_in, k, penguin_list, a):
		"""
		Arguments:
			F_self			:=	multiplicative strength factor of penguin self-propulsion
			F_in			:=	multiplicative strength factor of the penguin boundary force
			k				:=	spring constant relating to the repulsive force of overlapping penguins
			penguin_list	:=	python list containing all penguins in the system
			a 				:=	average radius of all penguins in the system
		"""

		critical_radius = 1.30 * a

		F_selfPropulsion = F_self * self.alignment

		F_repulsion = 0 

		if (self.boundary == True):

			F_boundary = F_in * self.alignment

		else:

			F_boundary = 0

		for i in range(len(penguin_list)):

			if (penguin_list[i] != self):

				r = self.get_distance(penguin_list[i])
				r_mag = np.linalg.norm(r)

				if (r < critical_radius):

					F_repulsion += -k * r

		return F_selfPropulsion + F_boundary + F_repulsion

	def net_torque(self, T_in, T_n, T_align, penguin_list, a):
		"""
		Arguments:
			T_in 			:=	multiplicative strength factor of the boundary torque term
			T_n				:=	multiplicative strength factor of the random torque term
			T_align			:=	multiplicative strength factor of the alignment torque term
			penguin_list	:=	python list containing all penguins in the system
			a 				:=	average radius of all penguins in the system
		"""

		critical_radius = 1.3 * a

		T_alignment = np.zeros(3)

		if (self.boundary == True):

			distance_sum = np.zeros(3)

			for i in range(len(penguin_list)):

				r = self.get_distance(penguin_list[i])
				r_mag = np.linalg.norm(r)

				if (penguin_list[i] != self) && (penguin_list[i].boundary == True):

					distance_sum += self.get_distance(penguin_list[i])

				if (penguin_list[i] != self) && (r_mag < critical_radius)):

					T_alignment += T_align * (self.alignment - penguin_list[i].alignment)

			exterior_bisector = distance_sum / np.linalg.norm(distance_sum)

			delta_theta = self.alignment - exterior_bisector

			T_boundary = T_in * delta_theta

		else:

			T_boundary = 0

		eta = random.uniform(-1.0,1.0)

		T_noise = eta * T_n

		return T_boundary + T_noise + T_alignment









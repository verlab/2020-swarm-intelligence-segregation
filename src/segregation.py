#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    _________                                         __  .__                  #
#   /   _____/ ____   ___________   ____   _________ _/  |_|__| ____   ____     #
#   \_____  \_/ __ \ / ___\_  __ \_/ __ \ / ___\__  \\   __\  |/  _ \ /    \    #
#   /        \  ___// /_/  >  | \/\  ___// /_/  > __ \|  | |  (  <_> )   |  \   #
#  /_______  /\___  >___  /|__|    \___  >___  (____  /__| |__|\____/|___|  /   #
#          \/     \/_____/             \/_____/     \/                    \/    #
#  																				#
#################################################################################
# Based on Matlab code of Vinicius Graciano Santos 
# Author: Paulo Rezeck
#################################################################################
# Cluster Segregation: dAA < dAB
#  Example: dAA = 3.0 and dAB = 5.0
# 
#  Cluster Aggregation: dAA > dAB
#  Example: dAA = 5.0 and dAB = 3.0
# 
#  Radial Segregation.: dAA(1) < dAB < dAA(2) < ... < dAA(n)
#  Example: dAA = [3.0, 5.0] and dAB = 4.0
# 
#  For convexity, all values must be larger than sqrt(3)/9
#################################################################################
import numpy as np
import math
from metric import *
import matplotlib.animation as animation

np.seterr(divide='ignore', invalid='ignore')
class Segregation(object):
	"""docstring for Segregation"""
	def __init__(self, alpha=1.5, ROBOTS=15, GROUPS=3, WORLD=40, dt=0.08, dAA=[7], dAB=20, noise_sensor=0.05, noise_actuation=0.05, seed=None, radius=1000, display_mode=False, which_metric='', DEAD_ROBOTS=[]):
		self.alpha = alpha
		self.ROBOTS = ROBOTS
		self.GROUPS = GROUPS
		self.DEAD_ROBOTS = DEAD_ROBOTS
		self.DEAD_ROBOTS.sort()
		self.WORLD = WORLD
		self.dt = dt
		self.dAA = np.array(dAA)
		self.dAB = dAB
		self.noise_sensor = noise_sensor
		self.noise_actuation = noise_actuation
		self.seed = seed
		self.RADIUS = radius
		self.display_mode = display_mode
		self.which_metric = which_metric
		# validation step
		self.validation()
		# Initialization
		self.setup()
		
	def validation(self):
		if self.ROBOTS % self.GROUPS != 0:
			print ("ROBOTS must be a multiple of GROUPS\n")
			quit()

		if (len(self.dAA) > 1.0) and (len(self.dAA) != self.GROUPS):
			print ("length(dAA) must be equal to GROUPS\n")
			quit()

		if any(self.dAA <= math.sqrt(3.0)/9.0) and self.dAB <= math.sqrt(3.0)/9.0:
			print ("Collective potential function is not strictly convex!\n")
			print ("dAA and dAB must be larger than sqrt(3)/9\n")
			quit()

		# if (self.which_metric == ''):
		# 	if (len(self.dAA) > 1.0):
		# 		#print "Robots will segregate to a radial configuration"
		# 		self.which_metric = 'radial'
		# 	else:	
		# 		#print "Robots will segregate to a cluster configuration"
		# 		self.which_metric = 'cluster'

	def setup(self):
		x = np.array(range(1, self.ROBOTS+1))
		i, j = np.meshgrid(x, x)

		gpr = float(self.GROUPS)/float(self.ROBOTS)
		
		AA = (np.floor(gpr*(i-1.0)) == np.floor(gpr*(j-1.0)))*1.0
		AB = (np.floor(gpr*(i-1.0)) != np.floor(gpr*(j-1.0)))*1.0
		
		# Vectorization of dAA and dAB.
		if len(self.dAA) == 1.0:
			self.const = np.multiply(self.dAA, AA) + np.multiply(self.dAB, AB)
		else:
			self.const = np.kron(np.diag(self.dAA), np.ones((int(self.ROBOTS/self.GROUPS), int(self.ROBOTS/self.GROUPS)))) + np.multiply(self.dAB, AB)
		
		np.random.seed(self.seed)
		self.q = self.WORLD * (np.random.uniform(low=0.0, high=1.0, size=(self.ROBOTS, 2)) - 0.5) # position
		np.random.seed(None)

		self.v = np.zeros((self.ROBOTS, 2)) # velocity

		self.metric_data = []
		# Choice which metric should the robots use
		if self.which_metric == 'cluster':
			self.metric = ClusterMetric(self.GROUPS, self.ROBOTS)
			print ("Using Cluster Metric!\n")	
		elif self.which_metric == 'radial':
			self.metric = RadialMetric(self.GROUPS, self.ROBOTS, self.const)
			print ("Using Radial Metric!\n")	
		elif self.which_metric == 'average':
			self.metric = AverageMetric(self.GROUPS, self.ROBOTS)		
			print ("Using Average Metric!\n")	
		elif self.which_metric == 'ncluster':
			self.metric = NClusterMetric(self.GROUPS, self.ROBOTS, self.dAA)		
			self.metric2 = AverageMetric(self.GROUPS, self.ROBOTS)		
			print ("Using NCluster Metric!\n")	

	
		self.fig = plt.figure()
		self.fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
		self.ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False,
							xlim=(-self.WORLD, self.WORLD), ylim=(-self.WORLD, self.WORLD))
		self.ax.grid(color='gray', linestyle='-', linewidth=0.1)
		self.cmap = plt.get_cmap('hsv')
		self.colors = [self.cmap(i) for i in np.linspace(0, 500/self.GROUPS, 500)]

		self.handler = []
		for i in range(self.GROUPS):
			self.particles, = self.ax.plot([], [], 'o', color=self.colors[i], ms=5)
			start = int(math.floor((i) * self.ROBOTS/self.GROUPS))
			stop = int(math.floor((i+1) * self.ROBOTS/self.GROUPS))
			self.particles.set_data(self.q[start:stop, 0], self.q[start:stop, 1])
			self.handler.append(self.particles)
		if self.display_mode:
			plt.ion()
			plt.show()


	def update(self):

		# Relative position among all pairs [q(j:2) - q(i:2)].
		xij = np.subtract(np.repeat(self.q[:,0], self.ROBOTS).reshape(self.ROBOTS, self.ROBOTS), self.q[:,0])
		yij = np.subtract(np.repeat(self.q[:,1], self.ROBOTS).reshape(self.ROBOTS, self.ROBOTS), self.q[:,1])
	
		# Relative velocity among all pairs [v(j:2) - v(i:2)].
		vxij = np.subtract(np.repeat(self.v[:,0], self.ROBOTS).reshape(self.ROBOTS, self.ROBOTS), self.v[:,0])
		vyij = np.subtract(np.repeat(self.v[:,1], self.ROBOTS).reshape(self.ROBOTS, self.ROBOTS), self.v[:,1])
	
		# Relative distance among all pairs.
		dsqr = xij**2 + yij**2
		
		# Add noise to sensor
		if self.noise_sensor != 0.00:
			s = (dsqr) * (self.noise_sensor)/3.0 + np.finfo(float).eps
			# Setup Normal-RNG to operate with max "noise" percent error with an error of 99.7% 3s = e
			dsqr = np.random.normal(dsqr, s)

		dist = np.sqrt(dsqr)
		# dist[dist > self.WORLD] = self.WORLD
	
		# Control equation.
		dU = np.multiply(self.alpha, (dist - self.const + 1.0/dist - self.const/dsqr))
		dU[dist > 2.0*self.RADIUS] = 0.0

		# for r in self.DEAD_ROBOTS:
		# 	dU[r,:] = 0.0
		# 	dU[:,r] = 0.0
		
		ax = np.multiply(-dU, xij)/dist - vxij # damping
		ay = np.multiply(-dU, yij)/dist - vyij # damping
	
		# a(i, :) -> acceleration input for robot i.
		self.a = np.array([np.nansum(ax, axis=1), np.nansum(ay, axis=1)]).T

		# Add noise to sensor
		if self.noise_actuation != 0.00:
			s = abs(self.a) * (self.noise_actuation)/3.0 + np.finfo(float).eps
			# Setup Normal-RNG to operate with max "noise" percent error with an error of 99.7% 3s = e
			self.a = np.random.normal(self.a, s)
	
	 	# simple taylor expansion.
		self.q = self.q + self.v*self.dt + self.a*(0.5*self.dt**2)
		self.v_n = self.v + self.a*self.dt
		self.v_n[(self.v_n[:,0] - self.v[:,0]) > 0.1] = 0.1
		self.v_n[(self.v_n[:,0] - self.v[:,0]) < -0.1] = -0.1
		self.v_n[(self.v_n[:,1] - self.v[:,1]) > 0.1] = 0.1
		self.v_n[(self.v_n[:,1] - self.v[:,1]) < -0.1] = -0.1
		self.v = self.v_n

		if self.which_metric == 'cluster':
			score = self.metric.compute(self.q, self.DEAD_ROBOTS)
			self.metric_data.append(score)
		elif self.which_metric == '':
			score = 0
			self.metric_data.append(score)

		m = np.array(abs(self.v))
		return m[:,0].mean() < 0.001 # stop
	
	def get_score(self):
		score = self.metric.compute(self.q)
		score2 = self.metric2.compute(self.q)
		self.metric_data.append([score, score2[0], score2[1]])

	def animate(self, i):
		if self.display_mode:
			# Update data for drawing.
			for i in range(self.GROUPS):
				start = int(math.floor((i) * self.ROBOTS/self.GROUPS))
				stop = int(math.floor((i+1) * self.ROBOTS/self.GROUPS))
				self.handler[i].set_data(self.q[start:stop, 0], self.q[start:stop, 1])
			return tuple(self.handler)
		else:
			print("Display is set to false!")
			return None


	def display(self):
		self.ax.clear()
		self.ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False,
								xlim=(-self.WORLD, self.WORLD), ylim=(-self.WORLD, self.WORLD))
		self.ax.grid(color='gray', linestyle='-', linewidth=0.1)
		self.ax.set_xlim([-self.WORLD, self.WORLD])
		self.ax.set_ylim([-self.WORLD, self.WORLD])
		self.ax.tick_params(axis='both', which='major', labelsize=6)
		self.ax.set_xlabel("X (meters)", fontsize=7)
		self.ax.set_ylabel("Y (meters)", fontsize=7)
		
		#self.ax.title("")
		plt.tight_layout()
		plt.gcf().subplots_adjust(bottom=0.2)

		for i in range(self.GROUPS):
			start = int(math.floor((i) * self.ROBOTS/self.GROUPS))
			stop = int(math.floor((i+1) * self.ROBOTS/self.GROUPS))
			robots = []
			for k in range(start, stop):
				if not k in self.DEAD_ROBOTS:
					robots.append(self.q[k,:])
			robots = np.array(robots)
			if (self.RADIUS < 1000):
				self.ax.plot(robots[:, 0], robots[:, 1], 'o', color='black', ms=2*self.RADIUS*1.4, fillstyle='none', alpha=0.6,  markeredgecolor=self.colors[i], markeredgewidth=0.1) # fig.dpi/72.
		for i in range(self.GROUPS):
				start = int(math.floor((i) * self.ROBOTS/self.GROUPS))
				stop = int(math.floor((i+1) * self.ROBOTS/self.GROUPS))
				robots = []
				for k in range(start, stop):
					if not k in self.DEAD_ROBOTS:
						robots.append(self.q[k,:])
				robots = np.array(robots)
				self.ax.plot(robots[:, 0], robots[:, 1], 'o', markerfacecolor=self.colors[i], fillstyle='none', ms=5, markeredgecolor='black', markeredgewidth=0.1)
		if self.display_mode:		
			plt.draw()
			plt.pause(0.0001)


	def screenshot(self, filename):
		plt.savefig(filename, dpi=100)

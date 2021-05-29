#!/usr/bin/env python
# -*- coding: utf-8 -*-

from segregation import Segregation
from termcolor import colored
import progressbar
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy as sp
import argparse


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h


# Setup - Getting user arguments
parser = argparse.ArgumentParser(description='Segregation using potential differential.')
parser.add_argument('--robots', help='Total number of robots.', default=50, type=int)
parser.add_argument('--groups', help='Total number of group of robots.', default=5, type=int)
parser.add_argument('--iterations', help='Total number of iterations on the control.', default=200, type=int)
# Reduce the steps if the robots are too shaky
parser.add_argument('--behavior', help='Choice between cluster (default) and radial', default="cluster", type=str)
parser.add_argument('--steps', help='Time-step on the control.', default=0.08, type=int)
parser.add_argument('--world', help='Size of the enviroment in meters', default=30, type=int)
parser.add_argument('--alpha', help='Control gain. See the paper for more information.', default=1.0, type=float)
parser.add_argument('--dAA', help='Same-type robot interaction factor. See the paper for more information.', default=2.0, type=float)
parser.add_argument('--dAB', help='Differente-type robot interaction factor. See the paper for more information.', default=5.0, type=float)
parser.add_argument('--noise_sensor', help='Add gaussian noise (max 1.0) on sensor model.', default=0.0, type=float)
parser.add_argument('--noise_actuation', help='Add gaussian noise (max 1.0) on actuator model.', default=0.0, type=float)
parser.add_argument('--sensing_radius', help='Limit the sensing radius.', default=100000000.0, type=float)
parser.add_argument('--seed', help='Random seed.', default=0, type=int)
args = parser.parse_args()


bar = progressbar.ProgressBar(maxval=args.iterations,
                              widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])

print(colored("[Starting numeric simulation...]", 'green'))

# Radial segregation (check dAA)
if args.behavior == 'radial':
	s = Segregation(ROBOTS=args.robots, GROUPS=args.groups, WORLD=args.world, dt=args.steps, alpha=args.alpha,  noise_sensor=args.noise_sensor, noise_actuation=args.noise_sensor, dAA=np.linspace(args.dAA, args.groups*args.dAA, args.groups), dAB=args.dAB, seed=args.seed, radius=args.sensing_radius, display_mode=True, which_metric='')
else:
	# Cluster segregation
	s = Segregation(ROBOTS=args.robots, GROUPS=args.groups, WORLD=args.world, dt=args.steps, alpha=args.alpha,  noise_sensor=args.noise_sensor, noise_actuation=args.noise_sensor, dAA=np.array([args.dAA]), dAB=args.dAB, seed=args.seed, radius=args.sensing_radius, display_mode=True, which_metric='')

bar.start()
for j in range(args.iterations):
    bar.update(j)
    s.update()
    s.display()
	# s.screenshot("log/image_{:06d}.png".format(j))

# s.screenshot("radial.png")
plt.close('all')

print(colored("[Finish]", 'green'))

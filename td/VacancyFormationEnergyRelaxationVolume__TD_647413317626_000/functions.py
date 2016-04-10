# Commonly used functions
# Author: Junhao Li <streaver91@gmail.com>
# Version: 0.0.2

import sys
from scipy.optimize import fmin
from scipy.interpolate import interp1d
import numpy as np
import datetime

def fmin_jh(
	func,
	x0,
	args = (),
	xtol = 0.0001,
	ftol = 0.0001,
	maxfun = None,
	full_output = 0,
	disp = 1,
	retall = 0,
	callback = None
):
	def fshifted(*args):
		argsNew = (args[0] + x0,) + args[1:]
		return func(*argsNew)
	res = fmin(
		fshifted,
		np.zeros_like(x0),
		args = args,
		xtol = xtol,
		ftol = ftol,
		maxfun = maxfun,
		full_output = full_output,
		disp = disp,
		retall = retall,
		callback = callback
	)
	resNew = (res[0] + x0,) + res[1:]
	return resNew

def fmax_jh(func, x0, args = (), xtol = 0.0001, ftol = 0.0001):
	def fshifted(x):
		return -func(x + x0)
	res = fmin(
		fshifted,
		np.zeros_like(x0),
		args = args,
		xtol = xtol,
		ftol = ftol
	)
	res[0] += x0
	return res

def clock(msg):
	print datetime.datetime.now().time(), msg

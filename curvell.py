import math
from scipy.optimize import minimize
import numpy as np


def ll3(b, probs, conc, sigma = 1000,weibull_param=[2,1]):
	b0 = b[0]
	b1 = b[1]
	b3 = b[2]
	if (b3 <= 0 | b3 > 1.001): return(-1e10)
	xi = np.exp(b0+b1*conc)
	alpha = 1+xi
	l = np.log(alpha)
	if (min(alpha)-b3 <=0): return(-1e10)
	wk = weibull_param[1]
	wl = weibull_param[2]*(wk/(wk-1))**(1/wk)
	#ll = MVN prior + weibull prior
	ll = -(b0**2 + b1**2)/(2*sigma**2) - (b3/wl)**wk + (wk-1)*np.log(b3) + \
		sum(probs*np.log(alpha-b3)) + sum((1-probs)*np.log(b3)) - sum(l)
	return(-ll)

def ll3_grad(b, probs, conc,sigma = 1000,weibull_param=[2,1]):
	b0 = b[0]
	b1 = b[1]
	b3 = b[2]
	xi = np.exp(b0+b1*conc)
	alpha = 1+xi
	d = (alpha - b3)
	m = probs*xi / d
	l = xi/alpha
	wk = weibull_param[1]
	wl = weibull_param[2]*(wk/(wk-1))**(1/wk)
	g1 = -b0/sigma**2 +sum(m) - sum(l)
	g2 = -b1/sigma**2 +sum(conc*m) - sum(conc*l)
	g4 = (wk-1)/b3 - (wk/b3)*(b3/wl)**wk - sum(probs/d) + sum((1-probs)/b3)
	return(np.array([-g1,-g2,-g4]))

def ll2(b, probs, conc, sigma = 1000):
	b0 = b[0]
	b1 = b[1]
	p_sum = sum(probs)
	p_conc_sum = sum(probs*conc)
	ll = -(b0**2 + b1**2)/(2*sigma**2) + b0*p_sum + b1*p_conc_sum - sum(np.log(1 + np.exp(b0+b1*conc)))
	return(-ll)

def ll2_grad(b, probs, conc, sigma = 1000):
	b0 = b[0]
	b1 = b[1]
	p_sum = sum(probs)
	p_conc_sum = sum(probs*conc)
	xi = np.exp(b0+b1*conc)
	l = xi/(xi+1)
	g1 = -b0/sigma**2 + sum(probs) - sum(l)
	g2 = -b1/sigma**2 + sum(conc*probs) - sum(conc*l)
	return(np.array([-g1,-g2]))

def fit_curve(probs, conc, curve_type = 'auto', method = 'BFGS'):
	#TODO: make auto method for choosing ll3 if low conc. mortality is above a threshold
	#TODO: create better guesses for b
	if curve_type in ["2", "ll2", 2]:
		b = np.array([1., 1.])
		res = minimize(ll2, b, args = (probs, conc), method = method, jac = ll2_grad )
	elif curve_type in ["3", "ll3", 3]:
		b = np.array([1., 1.])
		res = minimize(ll3, b, args = (probs, conc), method = method, jac = ll3_grad )
	return res



def MC_fit_curves(probs, conc, iter = 100, curve_type = 'auto', method = 'BFGS'):
	#TODO: 

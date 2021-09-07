import math
from scipy.optimize import minimize
import numpy as np


def ll3(b, probs, conc, sigma = 1000,weibull_param=[2,1]):
	b0 = b[0]
	b1 = b[1]
	b3 = b[2]
	if (b3 <= 1e-6 or b3 > 1 - 1e-6): return(1e10)
	xi = np.exp(b0+b1*conc)
	alpha = 1+xi
	l = np.log(alpha)
	if (min(alpha)-b3 <= 1e-6): return(1e10)
	wk = weibull_param[0]
	wl = weibull_param[1]*(wk/(wk-1))**(1/wk)
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
	wk = weibull_param[0]
	wl = weibull_param[1]*(wk/(wk-1))**(1/wk)
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
	#TODO: create better guesses for b
	#find lowest 20% of concentrations
	n_vals = round(0.2 * len(conc))
	idx = np.argpartition(conc, -n_vals)
	background_mort = sum(probs[idx[-n_vals:]])/n_vals
	if curve_type == 'auto':
		if background_mort > 0.15:
			curve_type = "ll3"
			med = (1. + background_mort)/2.
		else:
			curve_type = "ll2"
			med = 0.5
	high_idx = np.where(probs > med)[0]
	est_lc50 = conc[high_idx[-1]]
	if curve_type in ["2", "ll2", 2]:
		b = np.array([est_lc50, -1.])
		res = minimize(ll2, b, args = (probs, conc), method = method, jac = ll2_grad )
	elif curve_type in ["3", "ll3", 3]:
		b2 = np.array([est_lc50, -1., background_mort])
		res = minimize(ll3, b2, args = (probs, conc), method = method, jac = ll3_grad )
	return res


def MC_fit_curves(live_count, dead_count, conc, iter = 100, curve_type = 'auto', method = 'BFGS'):
	#TODO: 

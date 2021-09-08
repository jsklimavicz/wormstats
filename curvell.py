import math
from scipy.optimize import minimize
from scipy.special import gamma
import numpy as np
from corr_beta import corr_beta_vars

import matplotlib.pyplot as plt


def ll3(b, probs, conc, sigma = 1000,weibull_param=[2,1]):
	b0 = b[0]
	b1 = b[1]
	b3 = b[2]
	if (b3 <= 1e-10 or b3 > 1. ): return(1e10)
	xi = np.exp(b0+b1*conc)
	alpha = 1+xi
	l = np.log(alpha)
	if (min(alpha)-b3 <= 1e-10): return(1e10)
	wk = weibull_param[0]
	wl = weibull_param[1]*(wk/(wk-1))**(1/wk)
	#ll = MVN prior + weibull prior
	# ll = -(b0**2 + b1**2)/(2*sigma**2) - (b3/wl)**wk + (wk-1)*np.log(b3) + np.log(wk) + wk*np.log(wl) +\
	# 	sum(probs*np.log(alpha-b3)) + np.log(b3)*sum((1-probs)) - sum(l)
	ll = -(b0**2 + b1**2)/(2*(sigma**2)) + sum(probs*np.log(alpha-b3)) + np.log(b3)*sum((1-probs)) - sum(l)
	if b3 < 1 - 1e-9: 
		ll += - ((b3/wl)**wk) + (wk-1)*np.log(b3) #+ np.log(wk) - wk*np.log(wl)
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
	g1 = -b0/(sigma**2) + sum(m) - sum(l)
	g2 = -b1/(sigma**2) + sum(conc*m) - sum(conc*l)
	g4 = (wk-1)/b3 - (wk/b3)*((b3/wl)**wk) - sum(probs/d) + sum((1-probs)/b3)
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

def loglogit2(b, conc):
	return 1./(1.+ np.exp(-b[0] - conc * b[1]))

def loglogit3(b, conc):
	return 1 - (b[2]/(1.+ np.exp(-b[0] - conc * b[1])) + 1-b[2])

def fit_curve(probs, conc, curve_type = 'auto', method = 'BFGS'):
	#TODO: create better guesses for b1
	default_b1 = -1.
	n_vals = round(0.2 * len(conc))
	idx = np.argpartition(conc, n_vals)
	background_mort = sum(probs[idx[:n_vals]])/n_vals
	med = background_mort/2.
	if curve_type == 'auto':
		if background_mort > 0.15:
			curve_type = "ll3"
		else:
			curve_type = "ll2"
			med = 0.5
	high_idx = np.where(probs > med)[0]
	est_lc50 = conc[high_idx[-1]]
	if curve_type in ["2", "ll2", 2]:
		b2 = np.array([est_lc50, default_b1])
		return minimize(ll2, b2, args = (probs, conc), method = method, jac = ll2_grad)
	elif curve_type in ["3", "ll3", 3]:
		b3 = np.array([est_lc50, default_b1, background_mort])
		res = minimize(ll3, b3, args = (probs, conc), method = method, jac = ll3_grad)
		if not res.success:
			res = minimize(ll3, b3, args = (probs, conc), method = 'Nelder-Mead')
		return res
	elif curve_type.lower() in ["best", "aic"]:
		b2 = np.array([est_lc50, default_b1])
		res2 = minimize(ll2, b2, args = (probs, conc), method = method, jac = ll2_grad)
		b3 = np.array([res2.x[0], res2.x[1], background_mort])
		res3 = minimize(ll3, b3, args = (probs, conc), method = method, jac = ll3_grad)
		if not res3.success:
			res3 = minimize(ll3, b3, args = (probs, conc), method = 'Nelder-Mead')
		AIC2 = 4 - 2*res2.fun
		AIC3 = 6 - 2*res3.fun
		# print("AIC 2:", AIC2, "AIC3:", AIC3)
		print("AIC 2:", res2.fun, "AIC3:", res3.fun)
		return res2 if AIC2 <= AIC3 else res3


def MC_fit_curves(plate_ids, live_count, dead_count, conc, bs_iter = 100, curve_type = 'best', method = 'BFGS', scale = 0.5, rho = 0.25):
	#get correlated beta variables
	beta_probs = corr_beta_vars(plate_ids, live_count, dead_count, size = bs_iter, scale = scale, rho = rho)
	#possibility of having 2 or 3 parameters
	params = np.zeros((bs_iter, 3))
	iter_count = 0
	x = np.linspace(min(conc), max(conc), 101)
	# plt.plot(conc, live_count/(live_count + dead_count), 'g.')
	for row in beta_probs:
		# plt.plot(conc, row, '.')
		res = fit_curve(row, conc, curve_type = curve_type, method = method)
		# print(res)
		if len(res.x) == 3:
			params[iter_count] = res.x
		else:
			params[iter_count,0:2] = res.x
			params[iter_count,2] = 1
		y = loglogit3(params[iter_count], x)
		plt.plot(x, y)
		iter_count += 1
	# plt.show()
	return params
	



#plate_ids=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
#conc = np.array([-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7])
#live_count = np.array([11, 11, 12, 3, 3, 2, 0, 0, 0, 0, 9, 8, 7, 7, 5, 2, 1, 1, 0, 0, 11, 12, 14, 3, 5, 1, 0, 1, 0, 0, 9, 11, 5, 9, 3, 1, 0, 0, 0, 0])
#dead_count = np.array([5, 5, 8, 17, 14, 15, 14, 15, 18, 15, 5, 6, 5, 9, 10, 14, 18, 17, 19, 13, 4, 4, 6, 12, 10, 14, 13, 12, 19, 19, 4, 3, 10, 6, 11, 18, 14, 15, 16, 11])
#probs = dead_count/(dead_count + live_count)
#b = np.array([-1.2646771, 1.0533963, 0.7197049])
#b = np.array([-1.5646771, 1.1533963, 0.8197049])

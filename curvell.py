import math
from scipy.optimize import minimize
from scipy.special import gamma
import numpy as np
from corr_beta import corr_beta_vars
from joblib import load, Parallel, delayed
import matplotlib.pyplot as plt

class CI_finder:
	default_options = {"bs_iter": 1000, 
						"curve_type": 'll3', 
						"method": 'BFGS', 
						"scale": 0.5, 
						"rho": 0.25,
						"n_points": 101}
	def __init__(self, plate_ids, live_count, dead_count, conc):
		self.plate_ids = plate_ids
		self.live_count = live_count
		self.dead_count = dead_count
		self.conc = conc
		self.options = self.default_options
	def ll3(self, b, probs, conc, sigma = 1000,weibull_param=[2,1]):
		b0, b1, b3 = b
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
	def ll3_grad(self, b, probs, conc,sigma = 1000,weibull_param=[2,1]):
		b0, b1, b3 = b
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
	def ll2(self, b, probs, conc, sigma = 1000):
		b0, b1 = b
		p_sum = sum(probs)
		p_conc_sum = sum(probs*conc)
		ll = -(b0**2 + b1**2)/(2*sigma**2) + b0*p_sum + b1*p_conc_sum - sum(np.log(1 + np.exp(b0+b1*conc)))
		return(-ll)
	def ll2_grad(self, b, probs, conc, sigma = 1000):
		b0, b1 = b
		p_sum = sum(probs)
		p_conc_sum = sum(probs*conc)
		xi = np.exp(b0+b1*conc)
		l = xi/(xi+1)
		g1 = -b0/sigma**2 + sum(probs) - sum(l)
		g2 = -b1/sigma**2 + sum(conc*probs) - sum(conc*l)
		return(np.array([-g1,-g2]))
	def loglogit2(self, b, conc):
		return 1./(1.+ np.exp(-b[0] - conc * b[1]))
	def loglogit3(self, b, conc):
		return b[2]/(1.+ np.exp(b[0] + conc * b[1]))
	def ll3_find_EC(params, quantiles):
		quant2 = np.tile(np.array(quantiles), (len(params),1))
		params = np.reshape(np.repeat(np.array(ci.params[:,0:2]), len(quantiles)), 
			(ci.params[:,0:2].shape[0], ci.params[:,0:2].shape[1], len(quantiles)))
		return (np.log(1./quant2 - 1.)-params[:,0])/params[:,1]
	def fit_curve(self, probs):
		#TODO: better methods of b estimation
		#TODO: handle cases where EC50 is not going to be in the data :( 
		default_b1 = -1.
		n_vals = round(0.2 * len(self.conc))
		idx = np.argpartition(self.conc, n_vals)
		background_mort = sum(probs[idx[:n_vals]])/n_vals
		med = background_mort/2.
		if self.options["curve_type"].lower() == 'auto':
			if background_mort > 0.15:
				self.options["curve_type"] = "ll3"
			else:
				self.options["curve_type"] = "ll2"
				med = 0.5
		high_idx = np.where(probs > med)[0]
		est_lc50 = self.conc[high_idx[-1]]
		if self.options["curve_type"].lower() in ["2", "ll2", 2]:
			b2 = np.array([est_lc50, default_b1])
			return minimize(ll2, b2, args = (probs, self.conc), method = self.options["method"], jac = ll2_grad)
		elif self.options["curve_type"].lower() in ["3", "ll3", 3]:
			b3 = np.array([est_lc50, default_b1, background_mort])
			res = minimize(ll3, b3, args = (probs, self.conc), method = self.options["method"], jac = ll3_grad)
			if not res.success:
				res = minimize(ll3, b3, args = (probs, self.conc), method = 'Nelder-Mead')
			return res
		elif self.options["curve_type"].lower() in ["best", "aic"]:
			b2 = np.array([est_lc50, default_b1])
			res2 = minimize(ll2, b2, args = (probs, self.conc), method = self.options["method"], jac = ll2_grad)
			b3 = np.array([res2.x[0], res2.x[1], background_mort])
			res3 = minimize(ll3, b3, args = (probs, self.conc), method = self.options["method"], jac = ll3_grad)
			if not res3.success:
				res3 = minimize(ll3, b3, args = (probs, self.conc), method = 'Nelder-Mead')
			AIC2 = 4 - 2*res2.fun
			AIC3 = 6 - 2*res3.fun
			return res2 if AIC2 <= AIC3 else res3
	def bootstrap_CIs(self):
		#get correlated beta variables
		beta_probs = corr_beta_vars(self.plate_ids, 
									self.live_count, 
									self.dead_count, 
									size = self.options["bs_iter"], 
									scale = self.options["scale"], 
									rho = self.options["rho"])
		#possibility of having 2 or 3 parameters
		self.params = np.zeros((self.options["bs_iter"], 3))
		dict_list = Parallel(n_jobs=-1)(delayed(self.fit_curve)(row) for row in beta_probs)
		for iter_count, res in enumerate(dict_list):
			if len(res.x) == 3:
				self.params[iter_count] = res.x
			else:
				self.params[iter_count,0:2] = res.x
				self.params[iter_count,2] = 1.
	def calculate_curves(self):
		self.min_conc, self.max_conc = min(self.conc)-1, max(self.conc)+1
		self.x = np.linspace(self.min_conc, self.max_conc, self.options["n_points"])
		self.points = np.zeros((len(self.params), self.options["n_points"]))
		for iter_count, row in enumerate(self.params): self.points[iter_count] = loglogit3(row, self.x)
	def get_plot_CIs(self, quantiles = [.025, 0.5, 0.975], options=None):
		if options: self.update_options(options)
		self.quantiles = quantiles
		self.bootstrap_CIs()
		self.calculate_curves()
		self.quant = np.quantile(self.points, self.quantiles, interpolation='linear', axis = 0)
	def get_EC_CIs(self, EC = np.array([0.1, 0.5, 0.9]), CI_val=0.95, options=None):
		EC = 1-EC
		EC_vals = ll3_find_EC(self.params, EC)
		alpha = 1. - CI_val
		quantiles = np.array([alpha/2, 0.5, 1- alpha/2])
		EC_val_summary = np.quantile(EC_vals, quantiles, interpolation='linear', axis = 0)
		return np.transpose(EC_val_summary), 1 - EC
	def plot_CIs(self, low = .025, line = 0.5, high =  0.975, options=None):
		self.get_plot_CIs(quantiles = [low, line, high], options=options)
		plt.fill_between(self.x, self.quant[0], self.quant[2], alpha = 0.5)
		plt.plot(self.x, self.quant[1], 'r-')
		plt.plot(self.conc, self.live_count/(self.dead_count + self.live_count), 'g.')
		plt.show()
	def update_options(self, options):
		for k, v in options.items(): self.options[k] = v

	
# ci = CI_finder(plate_ids, live_count, dead_count, conc)
# ci.plot_CIs()
# ci.get_EC_CIs()


# plate_ids=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# conc = np.array([-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7])
# live_count = np.array([11, 11, 12, 3, 3, 2, 0, 0, 0, 0, 9, 8, 7, 7, 5, 2, 1, 1, 0, 0, 11, 12, 14, 3, 5, 1, 0, 1, 0, 0, 9, 11, 5, 9, 3, 1, 0, 0, 0, 0])
# dead_count = np.array([5, 5, 8, 17, 14, 15, 14, 15, 18, 15, 5, 6, 5, 9, 10, 14, 18, 17, 19, 13, 4, 4, 6, 12, 10, 14, 13, 12, 19, 19, 4, 3, 10, 6, 11, 18, 14, 15, 16, 11])
# probs = dead_count/(dead_count + live_count)
# b = np.array([-1.2646771, 1.0533963, 0.7197049])
# b = np.array([-1.5646771, 1.1533963, 0.8197049])

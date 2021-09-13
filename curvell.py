import math
from scipy.optimize import minimize
from scipy.special import gamma
import numpy as np
from corr_beta import corr_beta_vars
from joblib import load, Parallel, delayed
import matplotlib.pyplot as plt
import warnings
from statistics import median
from sklearn.neighbors import KernelDensity

def default_params_dict():
	param_dict = {"bs_iter": 500, 
						# "curve_type": 'll3', 
						"curve_type": 'best', 
						"method": 'BFGS', 
						"scale": 0.1, 
						"rho": 0.25,
						"n_points": 101,
						"CI_method": "HPDR"}
	return param_dict

class CI_finder:
	default_options = default_params_dict()

	def __init__(self, plate_ids, live_count, dead_count, conc, **kwargs):
		if plate_ids is None or live_count is None or dead_count is None or conc is None:
			raise ValueError("CI_finder requries plate_ids, live_count, dead_count, and conc values.")
		elif not (len(plate_ids)== len(live_count) and 
			len(plate_ids)==len(dead_count) and 
			len(plate_ids)==len(conc)):
			raise ValueError("CI_finder requries plate_ids, live_count, dead_count, and conc to all be the same length.")
		self.plate_ids = plate_ids
		self.live_count = live_count
		self.dead_count = dead_count
		self.conc = conc
		self.options = self.default_options
		#values created by calculation
		self.params = None
		self.points = None
		self.plot_quant = None

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

	def loglogit2(self, b, conc): return 1./(1.+ np.exp(-b[0] - conc * b[1]))

	def loglogit3(self, b, conc): return b[2]/(1.+ np.exp(b[0] + conc * b[1]))

	def ll3_find_EC(self, quantiles):
		quant2 = np.tile(np.array(quantiles), (len(self.params),1))
		params = np.reshape(np.repeat(np.array(self.params[:,0:2]), len(quantiles)), 
			(self.params[:,0:2].shape[0], self.params[:,0:2].shape[1], len(quantiles)))
		return (np.log(1./quant2 - 1.)-params[:,0])/params[:,1]

	def fit_curve(self, probs):
		#TODO: better methods of b estimation
		#TODO: handle cases where EC50 is not going to be in the data :( 
		def fit_ll3(b3, probs):
			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				res = minimize(self.ll3, b3, args = (probs, self.conc), method = self.options["method"], jac = self.ll3_grad)
			if not res.success:
				res = minimize(self.ll3, b3, args = (probs, self.conc), method = 'Nelder-Mead')
			return res
		default_b1 = median(self.conc)
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
			return minimize(self.ll2, b2, args = (probs, self.conc), method = self.options["method"], jac = self.ll2_grad)
		elif self.options["curve_type"].lower() in ["3", "ll3", 3]:
			b3 = np.array([est_lc50, default_b1, background_mort])
			res = fit_ll3(b3, probs)
			return res
		elif self.options["curve_type"].lower() in ["best", "aic"]:
			b2 = np.array([est_lc50, default_b1])
			res2 = minimize(self.ll2, b2, args = (probs, self.conc), method = self.options["method"], jac = self.ll2_grad)
			b3 = np.array([res2.x[0], res2.x[1], background_mort])
			res3 = fit_ll3(b3, probs)
			AIC2 = 4 - 2*res2.fun
			AIC3 = 6 - 2*res3.fun
			return res2 if AIC2 <= AIC3 else res3

	def bootstrap_CIs(self):
		if self.params is not None: return
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
		if self.points is not None: return
		self.min_conc, self.max_conc = min(self.conc)-1, max(self.conc)+1
		self.x = np.linspace(self.min_conc, self.max_conc, self.options["n_points"])
		self.points = np.zeros((len(self.params), self.options["n_points"]))
		for iter_count, row in enumerate(self.params): self.points[iter_count] = self.loglogit3(row, self.x)
	
	def get_plot_CIs(self, quantiles = [.025, 0.5, 0.975], options=None):
		if self.plot_quant is None: 
			if options: self.update_options(options)
			self.bootstrap_CIs()
			self.calculate_curves()
		self.plot_quant = np.quantile(self.points, quantiles, interpolation='linear', axis = 0)
	
	def get_CIs(self, EC = np.array([0.1, 0.5, 0.9]), CI_val=0.95, CI_method= "HPDR", *args, **kwargs):
		if self.params is None: self.bootstrap_CIs()
		EC = 1-EC
		EC_vals = self.ll3_find_EC(EC)
		alpha = 1. - CI_val
		return self.get_HPDR_CIs(EC_vals, alpha, *args, **kwargs) if CI_method == "HPDR" else self.get_EC_CIs(EC_vals, alpha, *args, **kwargs) 


	def get_EC_CIs(self, EC_vals, alpha, *args, **kwargs):
		quantiles = np.array([alpha/2, 0.5, 1- alpha/2])
		EC_val_summary = np.quantile(EC_vals, quantiles, interpolation='linear', axis = 0)
		return np.transpose(EC_val_summary) if len(EC_val_summary.shape)>1 else EC_val_summary

	def EC_kernel(self, EC_val_list):
		return KernelDensity(kernel ='gaussian', bandwidth = 0.5).fit(EC_val_list)

	def errfn(self, p, alpha, kde):
		from scipy import integrate
		def fn(x):
			x = np.array(x).reshape(1, -1)
			pdf = np.exp(kde.score_samples(x))
			return pdf if pdf > p else 0

		lb = min(self.conc) - 10
		ub = max(self.conc) + 10
		prob = integrate.quad(fn, lb, ub, limit=200)[0]
		return (prob + alpha - 1.)**2

	def get_HPDR_CIs(self, EC_vals, alpha, *args, **kwargs):
		from scipy.optimize import fmin
		for val_iter in EC_vals.T:
			val = np.reshape(val_iter, (len(val_iter),1))
			kde = self.EC_kernel(val)
			p = fmin(self.errfn, x0=0, args = (alpha, kde))
		exit()

	def get_HPDR_CIs(self, EC_vals, alpha, *args, **kwargs):
		from scipy.optimize import fmin
		for val_iter in EC_vals.T:
			val = np.reshape(val_iter, (len(val_iter),1))
			kde = self.EC_kernel(val)
			p = fmin(self.errfn, x0=0, args = (alpha, kde))
		exit()





	def get_EC50_CI(self, CI_val=0.95): return self.get_param_CI(0, CI_val)
	def get_slope_CI(self, CI_val=0.95): return self.get_param_CI(1, CI_val)
	def get_baseline_mort_CI(self, CI_val=0.95): return self.get_param_CI(2, CI_val)

	def get_param_CI(self, parameter, CI_val):
		if self.params is None: self.bootstrap_CIs()
		alpha = 1. - CI_val
		quantiles = np.array([alpha/2, 1- alpha/2])
		return np.quantile(self.params[:,parameter], quantiles, interpolation='linear', axis = 0)

	def plot_CIs(self, low = .025, line = 0.5, high =  0.975, options=None):
		self.get_plot_CIs(quantiles = [low, line, high], options=options)
		plt.fill_between(self.x, self.plot_quant[0], self.plot_quant[2], alpha = 0.5)
		plt.plot(self.x, self.plot_quant[1], 'r-')
		plt.plot(self.conc, self.live_count/(self.dead_count + self.live_count), 'g.')
		plt.show()
		return plt
	
	def update_options(self, options): 
		for k, v in options.items(): self.options[k] = v
	
	def reset_curves(self):
		'''
		This function clears out the params, points, and CI lines to allow one do new bootstrapping.
		'''
		self.params = None
		self.points = None
		self.plot_quant = None




def errfn(p, alpha, kde):
	from scipy import integrate
	def fn(x):
		x = np.array(x).reshape(1, -1)
		pdf = np.exp(kde.score_samples(x))
		return pdf if pdf > p else 0
	def get_bounds():
		def fn(x):
			x = np.array(x).reshape(1, -1)
			pdf = np.exp(kde.score_samples(x))
			return pdf - p
		lb, ub = -10, 10
		lb = fmin(fn, x0=lb)
		ub = fmin(fn, x0=ub)
		return lb, ub
	lb, ub = get_bounds()
	prob = integrate.quad(fn, lb, ub, limit=200)[0]
	return (prob + alpha - 1.)**2


def get_HPDR_CIs(kde, alpha, *args, **kwargs):
	return fmin(errfn, x0=0, args = (alpha, kde))

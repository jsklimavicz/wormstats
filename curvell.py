import math
from scipy.optimize import minimize
from scipy.special import gamma
from scipy.interpolate import CubicSpline
import numpy as np
from corr_beta import corr_beta_vars
import matplotlib.pyplot as plt
import warnings
from statistics import median, mean
from sklearn.neighbors import KernelDensity
from merlin_grapher import MerlinGrapher
from time import time
import utils

def default_params_dict():
	param_dict = {"BOOTSTRAP_ITERS": 1000,  
						"CURVE_TYPE": 'best', 
						"FIT_METHOD": 'BFGS', 
						"BETA_PRIOR": 0.1, 
						"RHO": 0.1,
						"N_POINTS": 151
				 }
	return param_dict

class CI_finder:	
	'''
	Performs the curve-fitting and associated math for dose-response analysis for the bioassay.
	'''
	default_options = default_params_dict()

	def __init__(self, plate_ids, live_count, dead_count, conc, options = None, **kwargs):
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
		self.n_trials = kwargs['n_trials']
		self.options = self.default_options
		if options:
			for k, v in options.items(): self.options[k] = v
		#values created by calculation
		self.params = None
		self.points = None
		self.plot_quant = None

	def ll3(self, b, probs, conc, sigma = 1000,weibull_param=[2,1]):
		'''
		Log-likelihood function of the three parameter dose-response curve 
							    b3
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein priors are b0, b1 ~ MVN(0, sigma*I2) and b2 ~ Weibull(weibull_param).
		'''
		b0, b1, b3 = b
		if (b3 <= 1e-10 or b3 > 1. ): return(1e10)
		xi = np.exp(b0+b1*conc)
		alpha = 1+xi
		l = np.log(alpha)
		if (min(alpha)-b3 <= 1e-10): return(1e10)
		wk = weibull_param[0]
		wl = weibull_param[1]*(wk/(wk-1))**(1/wk)
		ll = -(b0**2 + b1**2)/(2*(sigma**2)) + sum(probs*np.log(alpha-b3)) + np.log(b3)*sum((1-probs)) - sum(l)
		if b3 < 1 - 1e-9: 
			ll += - ((b3/wl)**wk) + (wk-1)*np.log(b3) #+ np.log(wk) - wk*np.log(wl)
		return(-ll)

	def ll3_grad(self, b, probs, conc,sigma = 1000,weibull_param=[2,1]):
		'''
		Gradient of the log-likelihood function of the three parameter dose-response curve 
							    b3
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein priors are b0, b1 ~ MVN(0, sigma*I2) and b2 ~ Weibull(weibull_param).
		'''
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
		'''
		Log-likelihood function of the two parameter dose-response curve 
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein prior is b0, b1 ~ MVN(0, sigma*I2).
		'''
		b0, b1 = b
		p_sum = sum(probs)
		p_conc_sum = sum(probs*conc)
		ll = -(b0**2 + b1**2)/(2*sigma**2) + b0*p_sum + b1*p_conc_sum - sum(np.log(1 + np.exp(b0+b1*conc)))
		return(-ll)

	def ll2_grad(self, b, probs, conc, sigma = 1000):
		'''
		Gradient of the log-likelihood function of the two parameter dose-response curve 
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein prior is b0, b1 ~ MVN(0, sigma*I2).
		'''
		b0, b1 = b
		p_sum = sum(probs)
		p_conc_sum = sum(probs*conc)
		xi = np.exp(b0+b1*conc)
		l = xi/(xi+1)
		g1 = -b0/sigma**2 + sum(probs) - sum(l)
		g2 = -b1/sigma**2 + sum(conc*probs) - sum(conc*l)
		return(np.array([-g1,-g2]))

	def loglogit2(self, b, conc): 
		'''
		Calculates the values of the two parameter dose-response curve 
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		given b = [b0, b1] and x.
		'''
		return 1./(1.+ np.exp(-b[0] - conc * b[1]))

	def loglogit3(self, b, conc): 
		'''
		Calculates the values of the three parameter dose-response curve 
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		given b = [b0, b1, b2] and x.
		'''
		return b[2]/(1.+ np.exp(b[0] + conc * b[1]))

	def ll3_find_LC(self, quantiles):
		'''
		Given the list self.params of parameters from bootstrapping and the list of quantiles, this 
		method returns the concentrations at which each of the quantiles is met. 
		'''
		quant2 = np.tile(np.array(quantiles), (len(self.params),1))
		params = np.reshape(np.repeat(np.array(self.params[:,0:2]), len(quantiles)), 
			(self.params[:,0:2].shape[0], self.params[:,0:2].shape[1], len(quantiles)))
		return (np.log(1./quant2 - 1.)-params[:,0])/params[:,1]

	#DECORATOR
	def update_options(func, *args, **kwargs): 
		'''
		Decorator to update the options for this class if the options keyword is included. 
		'''
		def inner(*args, **kwargs):
			if "options" in kwargs:
				for k, v in options.items(): self.options[k] = v
			func(*args, **kwargs)
		return inner

	def fit_curve(self, probs):
		'''
		Curve-fitting driver. 
		'''
		#TODO: better methods of b estimation
		#TODO: handle cases where EC50 is not going to be in the data :( 
		def fit_ll3(b3, probs):
			'''
			Inner function that attempts to fit the three-paremeter dose-response curve with the option-specified
			method for minimization. If this method fails, the more computationally-expensive but better-behaved
			Nelder-Mead method is used. 
			'''
			with warnings.catch_warnings(): #Ignores warnings for, e.g., including gradients for non-gradient methods. 
				warnings.simplefilter("ignore")
				res = minimize(self.ll3, b3, args = (probs, self.conc), method = self.options["FIT_METHOD"], jac = self.ll3_grad)
			if not res.success:
				res = minimize(self.ll3, b3, args = (probs, self.conc), method = 'Nelder-Mead')
			return res
		default_b1 = median(self.conc)
		n_vals = round(0.2 * len(self.conc))
		idx = np.argpartition(self.conc, n_vals)
		background_mort = sum(probs[idx[:n_vals]])/n_vals
		med = background_mort/2.
		if self.options["CURVE_TYPE"].lower() == 'auto':
			if background_mort > 0.15:
				self.options["CURVE_TYPE"] = "ll3"
			else:
				self.options["CURVE_TYPE"] = "ll2"
				med = 0.5
		high_idx = np.where(probs > med)[0]
		est_lc50 = self.conc[high_idx[-1]]
		if self.options["CURVE_TYPE"].lower() in ["2", "ll2", 2]:
			b2 = np.array([est_lc50, default_b1])
			return minimize(self.ll2, b2, args = (probs, self.conc), method = self.options["FIT_METHOD"], jac = self.ll2_grad)
		elif self.options["CURVE_TYPE"].lower() in ["3", "ll3", 3]:
			b3 = np.array([est_lc50, default_b1, background_mort])
			res = fit_ll3(b3, probs)
			return res
		elif self.options["CURVE_TYPE"].lower() in ["best", "aic"]:
			b2 = np.array([est_lc50, default_b1])
			res2 = minimize(self.ll2, b2, args = (probs, self.conc), method = self.options["FIT_METHOD"], jac = self.ll2_grad)
			b3 = np.array([res2.x[0], res2.x[1], background_mort])
			res3 = fit_ll3(b3, probs)
			AIC2 = 4 - 2*res2.fun
			AIC3 = 6 - 2*res3.fun
			return res2 if AIC2 <= AIC3 else res3

	def bootstrap_CIs(self):
		'''
		Driver for the bootstrapping process. Because the curvefitting is computationally
		costly but is embrassingly parallel, this is best done in parallel. 
		'''
		import multiprocessing
		from joblib import Parallel, delayed

		if self.params is not None: return
		#get correlated beta variables
		beta_probs = corr_beta_vars(self.plate_ids, 
									self.live_count, 
									self.dead_count, 
									size = self.options["BOOTSTRAP_ITERS"], 
									scale = self.options["BETA_PRIOR"], 
									null_scale = self.options["BETA_PRIOR_0"], 
									rho = self.options["RHO"])
		#calculate point data points while we have the beta parameters. 
		self.get_point_error_bars(beta_probs)
		#possibility of having 2 or 3 parameters
		self.params = np.zeros((self.options["BOOTSTRAP_ITERS"], 3))
		cpu_count = multiprocessing.cpu_count() 
		if self.options["BOOTSTRAP_ITERS"] > cpu_count : self.options["BOOTSTRAP_ITERS"] = cpu_count
		dict_list = Parallel(n_jobs=self.options["BOOTSTRAP_ITERS"])(delayed(self.fit_curve)(row) for row in beta_probs)
		for iter_count, res in enumerate(dict_list):
			if len(res.x) == 3:
				self.params[iter_count] = res.x
			else:
				self.params[iter_count,0:2] = res.x
				self.params[iter_count,2] = 1.
		
	def get_point_error_bars(self, beta_probs):
		'''
		Calculates the error bars for each data point based on the beta variables from the bootstrap samples. 
		'''
		lb = (1 - self.options["ERROR_BAR_CI"])/2.
		ub = 1-lb
		errors = np.quantile(beta_probs, [ub, lb], interpolation='linear', axis = 0)
		probs = self.dead_count / (self.dead_count+self.live_count)
		probs = np.tile(np.array(probs), (2,1))
		self.error_bars = np.abs(errors-probs)

	def calculate_curves(self):
		'''
		Calculates the value of the dose response curves at each value in self.x for each set of parameters. 
		'''
		if self.points is not None: return
		self.min_conc, self.max_conc = min(self.conc)-1, max(self.conc)+1
		lb, ub = math.floor(self.min_conc), math.ceil(self.max_conc)
		self.x = np.linspace(lb-1, ub+1, self.options["N_POINTS"])
		self.points = np.zeros((len(self.params), self.options["N_POINTS"]))
		for iter_count, row in enumerate(self.params): self.points[iter_count] = self.loglogit3(row, self.x)
		self.find_r2()
		# print(self.r2)

	def find_r2(self):
		'''
		Finds the r**2 value of the curve using the formula 
		r**2 = 1 - SS_{res}/SS_{tot}, where SS_{res} = sum(y_i - f_i)^2 for y_i being the response at
		concentration i and the f_i the value of the dose-response at i, and SS_{tot} = sum(y_i - y_bar)^2, 
		with y_bar being the average response across all concentrations. 
		'''
		# center_curve = np.quantile(self.points, 0.5, interpolation='linear', axis = 0)
		center_curve = np.median(self.points, axis = 0)
		spline = CubicSpline(self.x, center_curve)
		spline_vals = spline(self.conc)
		probs = self.live_count / (self.dead_count+self.live_count)
		start = 0
		end = 10
		sum_of_square_resids = sum((spline_vals - probs) ** 2)
		sum_of_square_nofit = sum((probs - mean(probs)) ** 2)
		self.r2 = 1- sum_of_square_resids/sum_of_square_nofit

	@update_options
	def get_plot_CIs(self, quantiles = [.025, 0.5, 0.975], options=None):
		'''
		Driver for calculating the confidence intervals for the plot.
		'''
		if self.plot_quant is None: 
			self.bootstrap_CIs()
			self.calculate_curves()
		self.plot_quant = utils.calc_ET_CI(self.points, CI_level = self.options["LC_CI"], resample = False)

	def get_CIs(self):
		'''
		Driver for calculating dose-response intervals. 
		'''
		if self.params is None: self.bootstrap_CIs()
		LC_VALUES = 1-self.options["LC_VALUES"]
		EC_vals = self.ll3_find_LC(LC_VALUES)
		EC_CI = self.get_EC_CIs(EC_vals, self.options["LC_CI"])
		# print(EC_CI)
		# exit()
		return np.power(2,EC_CI)

	def get_EC_CIs(self, EC_vals, CI_val):
		'''
		Finds the confidence intervals for doses. 
		'''
		func = utils.calc_ET_CI if self.options["CI_METHOD"].lower() in utils.ET_VARS else utils.calc_HPDI_CI
		EC_val_summary = func(EC_vals, CI_level = CI_val)
		# print(EC_val_summary)
		return EC_val_summary.T if EC_val_summary.ndim>1 else EC_val_summary

	def get_LC50_CI(self, CI_val=0.95): return  (- self.get_param_CI(0, CI_val)/self.get_param_CI(1, CI_val)).reshape((-1))
	def get_slope_CI(self, CI_val=0.95): return self.get_param_CI(1, CI_val).reshape((-1))
	def get_baseline_mort_CI(self, CI_val=0.95): return self.get_param_CI(2, CI_val).reshape((-1))

	def get_param_CI(self, parameter, CI_val):
		if self.params is None: self.bootstrap_CIs()
		# print(self.options["CI_METHOD"].lower())
		func = utils.calc_ET_CI if self.options["CI_METHOD"].lower() in utils.ET_VARS else utils.calc_HPDI_CI
		return func(self.params[:,parameter], CI_level = CI_val)

	def LC_kernel(self, LC_val = 0.5):
		'''
		Generates a KernelDensity object from an LC value.
		'''
		vals = self.ll3_find_LC(quantiles = LC_val)
		return KernelDensity(kernel ='gaussian', bandwidth = 0.5).fit(vals)

	def plot_CIs(self):
		lb = (1 - self.options["CURVE_CI"])/2
		ub = 1 - lb
		self.get_plot_CIs(quantiles = [lb, 0.5, ub])
		# print(self.options)
		merlin_plot = MerlinGrapher(x = self.x, 
									lb = self.plot_quant[0],
									ub = self.plot_quant[2],
									line = self.plot_quant[1],
									conc = self.conc,
									probs = self.live_count/(self.dead_count + self.live_count), 
									error_bars = self.error_bars,
									n_trials = self.n_trials,
									options = self.options)
		return merlin_plot

	def reset_curves(self):
		'''
		This function clears out the params, points, and CI lines to allow one do new bootstrapping.
		'''
		self.params = None
		self.points = None
		self.plot_quant = None



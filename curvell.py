import numpy as np
from curvell import CI_finder
import matplotlib.pyplot as plt


class Compound:
	def __init__(self, **kwargs):
		self.info = empty_cpmd_dict()
		for k, v in kwargs.items(): self.info[k] = v
		self.curve_data = None
	def fit_data(self, **kwargs):
		for k, v in kwargs.items(): 
			self.info[k] = v
		self.curve_data = CI_finder(**self.info)
	def get_CIs(self, EC = np.array([0.1, 0.5, 0.9]), CI_val=0.95, CI_method = "HPDR", options=None):
		if self.curve_data is None: self.fit_data()
		p = self.curve_data.get_CIs(EC = EC, CI_val=CI_val, CI_method = CI_method, options=options)
		return p
	def get_plot(self):
		ax = self.curve_data.plot_CIs()
		plt.show()

	def __iadd__(self, other):
		for k, v in self.info.items(): 
			if k == "name": pass
			elif k == "max_conc":
				self.info[k] = max(self.info[k], other.info[k])
			elif k == "n_trials":
				self.info[k] += 1
			elif k in ["conc", "live_count", "dead_count", "ctrl_mort"]:
				self.info[k] = numpy.concatenate(self.info[k], other.info[k])
			else:
				self.info[k] = [*self.info[k], *other.info[k]]
		self.curve_data = None


def empty_cpmd_dict():
	return {"name": None, #name of the compound
				"ids": None, #list of ID codes associated with the compound, of length number of data points.
				"test_dates": None, #list of testing dates associated with the compound, of length number of data points.
				"max_conc": None, #maximum concentration of all tests. 
				"n_trials": None, #Number of distinct complete reps. 
				"column_IDs": None, 
				"row_IDs": None, 
				"conc": None, #list of testing dates associated with the compound, of length number of data points.
				"live_count": None, #list of testing dates associated with the compound, of length number of data points.
				"dead_count": None,
				"plate_ids": None,
				"reps": None,
				"ctrl_mort": None,
				"unique_plate_ids": None}






# plate_ids=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# conc = np.array([-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7])
# live_count = np.array([11, 11, 12, 3, 3, 2, 0, 0, 0, 0, 9, 8, 7, 7, 5, 2, 1, 1, 0, 0, 11, 12, 14, 3, 5, 1, 0, 1, 0, 0, 9, 11, 5, 9, 3, 1, 0, 0, 0, 0])
# dead_count = np.array([5, 5, 8, 17, 14, 15, 14, 15, 18, 15, 5, 6, 5, 9, 10, 14, 18, 17, 19, 13, 4, 4, 6, 12, 10, 14, 13, 12, 19, 19, 4, 3, 10, 6, 11, 18, 14, 15, 16, 11])
# ch = Compound(name="malathion", id="1036-14A", conc=conc, live_count=live_count, plate_ids=plate_ids, dead_count=dead_count)
# ch.fit_data()
# ch.get_EC_CIs()

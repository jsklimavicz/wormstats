import numpy as np
from curvell import CI_finder
import matplotlib.pyplot as plt


class Compound:
	def __init__(self, **kwargs):
		self.info = {"name": None,
				"ids": None,
				"test_dates": None,
				"max_conc": None,
				"n_trials": None, 
				"conc": None,
				"live_count": None,
				"dead_count": None,
				"plate_ids": None}
		for k, v in kwargs.items(): self.info[k] = v
		self.curve_data = None
		self.data = None
	def fit_data(self, **kwargs):
		for k, v in kwargs.items(): self.info[k] = v
		self.curve_data = CI_finder(**self.info)
	def get_EC_CIs(self, EC = np.array([0.1, 0.5, 0.9]), CI_val=0.95, options=None):
		if self.curve_data is None: self.fit_data()
		return self.curve_data.get_EC_CIs(EC = EC, CI_val=CI_val, options=options)
	def get_plot(self):
		self.curve_data.plot_CIs()



# plate_ids=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# conc = np.array([-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7])
# live_count = np.array([11, 11, 12, 3, 3, 2, 0, 0, 0, 0, 9, 8, 7, 7, 5, 2, 1, 1, 0, 0, 11, 12, 14, 3, 5, 1, 0, 1, 0, 0, 9, 11, 5, 9, 3, 1, 0, 0, 0, 0])
# dead_count = np.array([5, 5, 8, 17, 14, 15, 14, 15, 18, 15, 5, 6, 5, 9, 10, 14, 18, 17, 19, 13, 4, 4, 6, 12, 10, 14, 13, 12, 19, 19, 4, 3, 10, 6, 11, 18, 14, 15, 16, 11])
# ch = Compound(name="malathion", id="1036-14A", conc=conc, live_count=live_count, plate_ids=plate_ids, dead_count=dead_count)
# ch.fit_data()
# ch.get_EC_CIs()

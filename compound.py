import numpy as np
from curvell import CI_finder
import matplotlib.pyplot as plt


class Compound:
	def __init__(self, **kwargs):
		self.info = {"name": None,
				"id": None,
				"max_conc": None,
				"n_trials": None, 
				"conc": None,
				"live_count": None,
				"dead_count": None,
				"plate_ids": None}
		for k, v in kwargs.items(): self.info[k] = v
	def fit_data(self, **kwargs):
		for k, v in kwargs.items(): self.info[k] = v
		self.curve_data = CI_finder(**self.info)
	def get_EC_CIs(self, EC = np.array([0.1, 0.5, 0.9]), CI_val=0.95, options=None):
		return self.curve_data.get_EC_CIs(EC = np.array([0.1, 0.5, 0.9]), CI_val=0.95, options=None)




ch = Compound(name="malathion", id="1036-14A", conc=conc, live_count=live_count, plate_ids=plate_ids, dead_count=dead_count)
ch.fit_data()
ch.get_EC_CIs()
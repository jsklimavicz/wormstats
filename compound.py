import numpy as np
from curvell import CI_finder
import matplotlib.pyplot as plt
from merlin_grapher import MerlinGrapher

class Compound:
	def __init__(self, *args, **kwargs):
		self.data = empty_cpmd_dict()
		for k, v in kwargs.items(): self.data[k] = v
		self.options = {}
		self.curve_data = None
		self.plot = None
	def fit_data(self, options):
		for k, v in options.items(): self.options[k] = v
		self.curve_data = CI_finder(**self.data, options = self.options)
	# def get_CIs(self, EC = np.array([0.1, 0.5, 0.9]), CI_val=0.95, CI_method = "HPDR", options=None):
	def get_CIs(self,  **kwargs):
		if self.curve_data is None: self.fit_data()
		# p = self.curve_data.get_CIs(EC = EC, CI_val=CI_val, CI_method = CI_method, options=options)
		p = self.curve_data.get_CIs(**kwargs)
		return p
	def make_plot(self):
		self.plot = self.curve_data.plot_CIs()
		self.plot.set_labels(title = self.data["name"])
		self.plot.plot_data()
		return self.plot.ax

	# def get_graph_ticks

	def __add__(self, other):
		other = self.remove_duplicates(other)
		if other is None: return Compound(**self.data)
		for k, v in self.data.items(): 
			if k == "name": pass
			elif k == "max_conc":
				self.data[k] = max(self.data[k], other.data[k])
			elif k == "n_trials":
				self.data[k] += 1
			elif k in ["conc", "live_count", "dead_count", "ctrl_mort"]:
				self.data[k] = np.array([*self.data[k], *other.data[k]])
			else:
				print(k)
				print(self.data[k], other.data[k])
				self.data[k] = [*self.data[k], *other.data[k]]
		return Compound(**self.data)

	def __sub__(self, uid_list):
		for k, v in self.data.items(): 
			if k == "name": pass
			elif k in ["conc", "live_count", "dead_count", "ctrl_mort", "column_IDs", "row_IDs", "plate_ids", "reps", "ctrl_mort", "unique_plate_ids"]:
				del self.data[k][uid_list]
		date, plate, row, ids = [list(d) for d in zip(*[[i for i in c.split('_')] for c in self.data["unique_plate_ids"]])]
		self.data["ID"] = list(set(ids))
		self.data["test_dates"] = list(set(date))
		self.data["n_trials"] = list(set(self.data["n_trials"]))
		self.data["max_conc"] = max(self.data["conc"])
		return Compound(**self.data)

	def remove_duplicates(self, other):
		dup_uniqueIDs = [uid for uid in other.data["unique_plate_ids"] if uid in self.data["unique_plate_ids"]]
		if len(dup_uniqueIDs) > 0: #then duplicates exist
			dup_loc = [i for i, x in enumerate(other.data["unique_plate_ids"]) if x in dup_uniqueIDs]
			if len(dup_loc) == len(self.data["unique_plate_ids"]): return None #all values are duplicates
			else: return other - dup_loc
		else: return other

	def test_print(self):
		print(self.data["name"])
		for k, v in self.data.items(): 
			print("      ", k, ": ", "(", type(v), ")", v, sep = "")
		if self.curve_data is not None: print("Curve data found.")
		if self.plot is not None: print("Plot found.")

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

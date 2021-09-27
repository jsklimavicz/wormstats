import numpy as np
from curvell import CI_finder
import matplotlib.pyplot as plt
from merlin_grapher import MerlinGrapher

class Compound:
	'''
	Class to hold all information for a single compound, from the name, to live/dead counts, 
	to the options for graphing, etc.
	'''
	def __init__(self, *args, **kwargs):
		'''
		Initiation of this class should ideally contain as much data as possible from the bioassay. 
		See the empty_cpmd_dict function in this module for potential keywords.
		'''
		self.data = empty_cpmd_dict()
		for k, v in kwargs.items(): self.data[k] = v
		self.data["max_conc"] = max(self.data["conc"])
		self.data["min_conc"] = min(self.data["conc"])
		self.options = {}
		self.curve_data = None
		self.plot = None
	def fit_data(self, options):
		'''
		Sets the curve data for this compounds using class CI_finder in module curvell
		'''
		for k, v in options.items(): self.options[k] = v
		# print(self.options)
		if self.curve_data is None: self.curve_data = CI_finder(**self.data, options = self.options)
	def get_LC_CIs(self):
		'''
		Calculates and returns confidence intervals for LC values. 
		'''
		return self.curve_data.get_CIs()
	def make_plot(self):
		'''
		Generates a graph of the dose-response data. The graph is initiated in the CI_finder class plot_CIs method, 
		as a MerlinGrapher class in the merlin_grapher module. 
		Returns the axis array. 
		'''
		self.plot = self.curve_data.plot_CIs()
		self.plot.set_labels(title = self.data["name"])
		self.plot.plot_data()
		# return self.plot.ax

	def reset_curves(self):
		self.curve_data = None
		CA = Compound(**self.data)
		CA.options = self.options
		return CA


	def __add__(self, other):
		'''
		Overloaded + operator to join two compound classes that have the same compound name. 
		'''
		other = self.remove_duplicates(other)
		if other is None: 
			cmpd = Compound(**self.data)
			cmpd.options = self.options
			cmpd.curve_data = self.curve_data
			return cmpd
		#If we get to this point, then we need to calculate new curves. 
		for k, v in self.data.items(): 
			if k == "name": pass
			elif k == "max_conc":
				self.data[k] = max(self.data[k], other.data[k])
			elif k == "min_conc":
				self.data[k] = min(self.data[k], other.data[k])
			elif k == "n_trials":
				self.data[k] += 1
			elif k in ["conc", "live_count", "dead_count", "ctrl_mort"]:
				self.data[k] = np.array([*self.data[k], *other.data[k]])
			else:
				self.data[k] = [*self.data[k], *other.data[k]]
		return Compound(**self.data)

	def __sub__(self, uid_list):
		'''
		Overloaded - operator to remove entries in the Compound class based on unique ID, to make sure, e.g., that 
		a specific trial is not included more than once. 
		'''

		# print(uid_list)
		# print(self.data)
		for k, v in self.data.items(): 
			if k == "name": pass
			elif k in ["conc", "live_count", "dead_count", "ctrl_mort"]:
				self.data[k] = np.delete(self.data[k], np.array(uid_list))
			elif k in ["column_IDs", "row_IDs", "plate_ids", "reps", "unique_plate_ids"]:
				delete_multiple_element(self.data[k], uid_list)
			# print(self.data[k])
		# print(self.data)
		date, plate, row, ids = [list(d) for d in zip(*[[i for i in c.split('_')] for c in self.data["unique_plate_ids"]])]
		self.data["ID"] = list(set(ids))
		self.data["test_dates"] = list(set(date))
		self.data["n_trials"] = len(list(set(self.data["test_dates"])))
		self.data["max_conc"] = max(self.data["conc"])
		self.data["min_conc"] = min(self.data["conc"])
		return Compound(**self.data)

	def remove_duplicates(self, other):
		'''
		Finds 
		'''
		dup_uniqueIDs = [uid for uid in other.data["unique_plate_ids"] if uid in self.data["unique_plate_ids"]]
		# print(len(self.data['column_IDs']))
		r = list(range(len(self.data['column_IDs'])))
		# print(len(other.data["unique_plate_ids"]))
		# print(len(dup_uniqueIDs))
		if len(dup_uniqueIDs) > 0: #then duplicates exist
			dup_loc = [i for i, x in enumerate(other.data["unique_plate_ids"]) if x in dup_uniqueIDs]
			# print(r)
			# print(dup_loc)
			r = delete_multiple_element(r, dup_loc)
			# print(r)
			if len(dup_loc) == len(other.data["unique_plate_ids"]): return None #all values are duplicates
			else: return other - dup_loc
		else: return other

	def test_print(self):
		print(self.data["name"])
		for k, v in self.data.items(): 
			print("      ", k, ": ", "(", type(v), ")", v, sep = "")
		if self.curve_data is not None: print("Curve data found.")
		if self.plot is not None: print("Plot found.")

	def saveable_cmpd(self): 
		'''
		Makes a lower memory version of the compound for saving purposes. Basically, we want to save the fitted dose-response
		curve parameters because these are expensive to claculate, but don't need to save the plot, or the fitted data because 
		these values are cheap to calculate but take up a lot of storage. 
		'''
		no_plt_cmpd = Compound(**self.data)
		no_plt_cmpd.options = self.options
		no_plt_cmpd.curve_data = self.curve_data
		#Remove the really storage-pricy stuff.
		no_plt_cmpd.curve_data.plot_quant = None
		no_plt_cmpd.curve_data.points = None
		return no_plt_cmpd

def empty_cpmd_dict():
	return {"name": None, #name of the compound
				"ids": None, #list of ID codes associated with the compound, of length number of data points.
				"test_dates": None, #list of testing dates associated with the compound, of length number of data points.
				"max_conc": None, #maximum concentration of all tests. 
				"min_conc": None,
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

def delete_multiple_element(list_object, indices):
    indices = sorted(indices, reverse=True)
    for idx in indices:
        if idx < len(list_object):
            list_object.pop(idx)
    return list_object

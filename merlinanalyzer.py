import numpy as np
import pandas as pd
from curvell import CI_finder
import matplotlib.pyplot as plt
from compound import Compound
import math

class MerlinAnalyzer:

	class MerlinData:
		def __init__(self, *args, **kwargs):
			self.params = default_params_dict()

	def __init__(self, *args, **kwargs):
		self.data = {}

	# def read_cmpd_key(self, filename):


	def read_new_data(self, filename):
		new_data = pd.read_csv(filename, header = 0)
		cmpd_data = new_data[new_data["ppm"] != 0]
		ctrl_data = new_data[new_data["ppm"] == 0]
		ctrl_live_sum = sum(np.array(ctrl_data["Live"].tolist()))
		ctrl_dead_sum = sum(np.array(ctrl_data["Dead"].tolist()))
		ctrl_ave_mort = (ctrl_dead_sum *1.)/(ctrl_dead_sum + ctrl_live_sum *1.)
		cmpd_data["ctrl_mort"] = ctrl_ave_mort
		self.new_data = cmpd_data

			
	def read_key(self, filename):
		this_data = pd.read_csv(filename, header = 0)

	def process_compounds(self):
		for cmpd_id in this_data["Ref ID"].unique():
			print(cmpd_id)
			cmpd_data = this_data[this_data["Ref ID"] == cmpd_id]
			cmpd = Compound(name = cmpd_id, 
							ids = cmpd_id, 
							test_dates = cmpd_data["Date"].unique(),
							max_conc = max(cmpd_data["ppm"].tolist()),
							n_trials = len(cmpd_data["Row"].unique()),
							column_IDs = cmpd_data["Column"].tolist(),
							row_IDs = cmpd_data["Row"].tolist(),
							conc = np.log(np.array(cmpd_data["ppm"].tolist()))/math.log(2),
							live_count = np.array(cmpd_data["Live"].tolist()),
							dead_count = np.array(cmpd_data["Dead"].tolist()),
							plate_ids = cmpd_data["Plate"].tolist(),
							reps = cmpd_data["Rep"].tolist())
			cmpd.fit_data()
			EC = cmpd.get_CIs(CI_method = "equal_tail")
			print(EC)
			# EC = cmpd.get_CIs(CI_method = "HPDR")
			# print(EC)
			cmpd.curve_data.plot_CIs()
		return

	def read_old_data(self, filename):
		return

	def save_json(self, json_filename):
		return

	def save_csv(Self, csv_filename):
		return

	def compare_LC(self, cmpd1, cmpd2, LCval = 50, CI = 0.95):
		return



if __name__ == "__main__":
	MA = MerlinAnalyzer()
	MA.read_new_data("./test.csv")

import numpy as np
import pandas as pd
from curvell import CI_finder
import matplotlib.pyplot as plt
from compound import Compound
import math
from copy import deepcopy

class MerlinAnalyzer:

	def __init__(self, *args, **kwargs):
		self.cmpd_data = {}

	# def read_cmpd_key(self, filename):


	def column_name_modifier(self, filename):
		new_data = pd.read_csv(filename, header = 0)
		column_names = list(new_data.columns)
		new_col_names = {}
		for col_name in column_names:
			if col_name.lower() in ["alive", 'live', 'living']: new_col_names[col_name] = "Live"
			elif col_name.lower() in ["dead"]: new_col_names[col_name] = "Dead"
			elif col_name.lower() in ["total", 'count', 'sum']: new_col_names[col_name] = "Count"
			elif col_name.lower() in ["ref.id", 'ref id', 'ref', 'id']: new_col_names[col_name] = "ID"
			elif 'col' in col_name.lower(): new_col_names[col_name] = "Column"
			elif 'row' in col_name.lower(): new_col_names[col_name] = "Row"
			elif col_name.lower() in ["ppm", "conc", "concentration"]: new_col_names[col_name] = "Conc"
			elif col_name.lower() in ["plate"]: new_col_names[col_name] = "Plate"
			elif col_name.lower() in ["rep"]: new_col_names[col_name] = "Rep"
			elif col_name.lower() in ["class"]: new_col_names[col_name] = "Class"
			elif col_name.lower() in ["date", 'day']: new_col_names[col_name] = "Date"
			else: new_col_names[col_name] = col_name
		new_data.rename(columns=new_col_names, inplace=True)
		return new_data

	def read_new_data(self, filename, *args, **kwargs):
		new_data = self.column_name_modifier(filename)
		cmpd_data = deepcopy(new_data[new_data["Conc"] != 0])
		ctrl_data = new_data[new_data["Conc"] == 0]
		ctrl_live_sum = sum(np.array(ctrl_data["Live"].tolist()))
		ctrl_dead_sum = sum(np.array(ctrl_data["Dead"].tolist()))
		ctrl_ave_mort = (ctrl_dead_sum *1.)/(ctrl_dead_sum + ctrl_live_sum *1.)
		cmpd_data["ctrl_mort"] = ctrl_ave_mort
		# cmpd_data.assign(ctrl_mort=ctrl_ave_mort)
		self.new_data = cmpd_data

	def read_merge_key(self, filename, *args, **kwargs):
		key = self.column_name_modifier(filename)[["Compound", "ID", "Class"]]
		self.new_data = self.new_data.merge(key, how='left', on='ID')

	def process_compounds(self, *args, **kwargs):
		for cmpd_id in self.new_data["Compound"].unique():
			print(cmpd_id)
			cmpd_data = self.new_data[self.new_data["Compound"] == cmpd_id].copy()
			print(cmpd_data)
			cmpd = Compound(name = cmpd_id, 
							ids = cmpd_id, 
							test_dates = cmpd_data["Date"].unique(),
							max_conc = max(cmpd_data["Conc"].tolist()),
							n_trials = len(cmpd_data["Row"].unique()),
							column_IDs = cmpd_data["Column"].tolist(),
							row_IDs = cmpd_data["Row"].tolist(),
							conc = np.log(np.array(cmpd_data["Conc"].tolist()))/math.log(2),
							live_count = np.array(cmpd_data["Live"].tolist()),
							dead_count = np.array(cmpd_data["Dead"].tolist()),
							plate_ids = cmpd_data["Plate"].tolist(),
							reps = cmpd_data["Rep"].tolist(),
							ctrl_mort = np.array(cmpd_data["ctrl_mort"].tolist()))
			cmpd.fit_data()
			EC = cmpd.get_CIs(CI_method = "equal_tail")
			print(EC)
			# EC = cmpd.get_CIs(CI_method = "HPDR")
			# print(EC)
			cmpd.get_plot()
		return

	def read_old_data(self, filename, *args, **kwargs):
		return

	def save_archive(self, filename, *args, **kwargs):
		return

	def save_csv(self, filename, *args, **kwargs):
		return

	def save_plots(self, image_dir, *args, **kwargs):
		return

	def compare_LC(self, cmpd1, cmpd2, LCval = 50, CI = 0.95, *args, **kwargs):
		return

	def full_process(self, 
						new_datafile = None, 
						old_datafile = None, 
						key_file = None,
						csv_outfile = None, 
						archive_outfile = None,
						image_dir = None, 
						*args, 
						**kwargs):
		if new_datafile: self.read_new_data(new_datafile, *args, **kwargs)
		if old_datafile: self.read_old_data(old_datafile, *args, **kwargs)
		if key_file: self.read_merge_key(key_file, *args, **kwargs)
		self.process_compounds(*args, **kwargs)

		if csv_outfile: self.save_csv(csv_outfile, *args, **kwargs)
		if archive_outfile: self.save_archive(archive_outfile, *args, **kwargs)
		if image_dir: self.save_plots(image_dir, *args, **kwargs)



if __name__ == "__main__":
	MA = MerlinAnalyzer()
	MA.full_process(new_datafile = "./test.csv", key_file = "./key.csv")
	# MA.read_new_data("./test.csv")
	# MA.process_compounds() 

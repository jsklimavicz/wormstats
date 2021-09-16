import numpy as np
import os
import pandas as pd
from curvell import CI_finder
import matplotlib.pyplot as plt
from compound import Compound
import math
from copy import deepcopy
import pickle
import hashlib
import hmac

class MerlinAnalyzer:

	archivefilename = "merlin_bioassay_archive_data.pickle"
	picklesha1hash = ".picklehash"

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

	def read_new_data(self, filename, key_file):
		new_data = self.column_name_modifier(filename)
		cmpd_data = deepcopy(new_data[new_data["Conc"] != 0])
		ctrl_data = new_data[new_data["Conc"] == 0]
		ctrl_live_sum = sum(np.array(ctrl_data["Live"].tolist()))
		ctrl_dead_sum = sum(np.array(ctrl_data["Dead"].tolist()))
		ctrl_ave_mort = (ctrl_dead_sum *1.)/(ctrl_dead_sum + ctrl_live_sum *1.)
		cmpd_data["ctrl_mort"] = ctrl_ave_mort
		if key_file is not None: key = self.read_merge_key(key_file)
		cmpd_data = cmpd_data.merge(key, how='left', on='ID')
		return self.process_compounds(cmpd_data)

	def read_merge_key(self, filename):
		return self.column_name_modifier(filename)[["Compound", "ID", "Class"]]
		
	def read_archive(self, filepath):
		recvd_digest, pickled_data = data.split(' ')
		new_digest = hmac.new(b'shared-key', pickled_data, hashlib.sha1).hexdigest()
		if recvd_digest != new_digest:
			print('Integrity check of pickled data has failed. Old data cannot be loaded.')
		else:
			unpickled_data = pickle.loads(pickled_data)

	def merge_old_new(self, new_datafile, old_datafile, key_file):
		if old_datafile is not None: self.cmpd_data = self.read_archive(old_datafile)
		if new_datafile is not None:
			new_cmpd_dict = self.read_new_data(new_datafile, key_file)

			#merge
			for k, v, in new_cmpd_dict.items():
				if k in self.cmpd_data: self.cmpd_data[k] += new_cmpd_dict[k]
				else: self.cmpd_data[k] = v

	def save_archive(self, filepath, *args, **kwargs):
		pickle_data = pickle.dumps(self.cmpd_data)
		digest =  hmac.new(b'shared-key', pickle_data, hashlib.sha1).hexdigest()
		header = '%s' % (digest)
		with open(os.path.join(filepath, self.picklesha1hash), 'w') as file:
			file.write(header)
		with open(os.path.join(filepath, self.archivefilename), 'wb') as file:
			file.write(pickle_data)
		
	def process_compounds(self, new_data, *args, **kwargs):
		new_compound_dict = {}
		for cmpd_id in new_data["Compound"].unique():
			cmpd_data = new_data[new_data["Compound"] == cmpd_id].copy()
			new_compound_dict[cmpd_id] = Compound(name = cmpd_id, 
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
							ctrl_mort = np.array(cmpd_data["ctrl_mort"].tolist()),
							*args, **kwargs)
		return new_compound_dict

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

		self.merge_old_new(new_datafile, old_datafile, key_file)
		# self.process_compounds(*args, **kwargs)

		for cmpd_name, cmpd in self.cmpd_data.items():
			print(cmpd.data["name"])
			cmpd.fit_data()
			EC = cmpd.get_CIs(CI_method = "equal_tail")
			print(EC)
		# # EC = cmpd.get_CIs(CI_method = "HPDR")
		# # print(EC)
			# ax = cmpd.plot_curve()
			# plt.show()

		if csv_outfile: self.save_csv(csv_outfile, *args, **kwargs)
		if archive_outfile: self.save_archive(archive_outfile, *args, **kwargs)
		if image_dir: self.save_plots(image_dir, *args, **kwargs)

if __name__ == "__main__":
	MA = MerlinAnalyzer()
	MA.full_process(new_datafile = "./test.csv", key_file = "./key.csv", archive_outfile=".")
	# MA.read_new_data("./test.csv")
	# MA.process_compounds() 

import numpy as np
import os
import sys
import pandas as pd
from curvell import CI_finder
import matplotlib.pyplot as plt
from compound import Compound
import math
from copy import deepcopy
import pickle
import hashlib
import hmac
from merlin_grapher import MerlinGrapher
import utils
from time import time

class MerlinAnalyzer:

	archivefilename = "merlin_bioassay_archive_data.pickle"
	picklesha1hash = ".picklehash"
	sha_key = b"merlin-data"

	def __init__(self, *args, **kwargs):
		self.cmpd_data = {}
		self.options = utils.parse_config_file()

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
		if key_file is not None: key = self.read_key(key_file)
		cmpd_data = cmpd_data.merge(key, how='left', on='ID')
		return self.process_compounds(cmpd_data)

	def read_key(self, filename):
		return self.column_name_modifier(filename)[["Compound", "ID", "Class"]]
		
	def read_archive(self, filepath):
		with open(os.path.join(filepath, self.picklesha1hash), 'r') as file:
			pickle_hash = file.read().strip()
		with open(os.path.join(filepath, self.archivefilename), 'rb') as file:
			pickled_data = file.read()
		digest =  hmac.new(self.sha_key, pickled_data, hashlib.sha1).hexdigest()
		if pickle_hash == digest:
			unpickled_data = pickle.loads(pickled_data)
			return unpickled_data
		else:
			print('Pickled data as been compromised. Old data cannot be loaded.')

	def merge_old_new(self, new_datafile, archive_path, key_file):
		if os.path.exists(os.path.join(archive_path,self.archivefilename)): 
			self.cmpd_data = self.read_archive(archive_path)

		if new_datafile is not None:
			new_cmpd_dict = self.read_new_data(new_datafile, key_file)
			#merge
			for k, v, in new_cmpd_dict.items():
				if k in self.cmpd_data: self.cmpd_data[k] = self.cmpd_data[k] + new_cmpd_dict[k]
				else: self.cmpd_data[k] = v
				# self.cmpd_data[k].test_print()

	def save_archive(self, filepath, *args, **kwargs):
		saveable_lib = {}
		for k, v in self.cmpd_data.items():
			saveable_lib[k] = self.cmpd_data[k].saveable_cmpd()
		pickle_data = pickle.dumps(saveable_lib)
		digest =  hmac.new(self.sha_key, pickle_data, hashlib.sha1).hexdigest()
		header = '%s' % (digest)
		with open(os.path.join(filepath, self.picklesha1hash), 'w') as file:
			file.write(header)
		with open(os.path.join(filepath, self.archivefilename), 'wb') as file:
			file.write(pickle_data)
		
	def process_compounds(self, new_data, *args, **kwargs):
		new_compound_dict = {}
		for cmpd_id in new_data["Compound"].unique():
			cmpd_data = new_data[new_data["Compound"] == cmpd_id].copy()
			unique_ids = ["_".join([x,str(y),str(z),w]) for x,y,z,w in zip(cmpd_data["Date"].tolist(), 
																		cmpd_data["Plate"].tolist(), 
																		cmpd_data["Row"].tolist(),
																		cmpd_data["ID"].tolist())]
			new_compound_dict[cmpd_id] = Compound(name = cmpd_id, 
							ids = cmpd_data["ID"].unique(), 
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
							unique_plate_ids = unique_ids,
							*args, **kwargs)
		return new_compound_dict

	def save_csv(self, filename, header, body):
		with open(filename, 'w') as file:
			file.write(",".join(header) + "\n")
			file.write("\n".join(body))

	def generate_csv_data_lines(self, header):
		def format_lc_val(val, sig = 3):
			if val < 1e-2 or val >1e3: 
				return f"{(round(val,sig)):.2e}"
			else: return f"{(round(val,sig)):.3g}"


		def CI_to_string(ci1, ci2, sig = 3):
			return '"[' + format_lc_val(ci1,sig) +", "+ format_lc_val(ci2,sig) + ']"'

		def format_LC_to_CSV(CI_array, sig = 3):
			output = []
			for line in CI_array:
				output.append(format_lc_val(line[1],sig))
				output.append(CI_to_string(line[0],line[2], sig = sig))
			return output

		output = []
		for cmpd_name, cmpd in self.cmpd_data.items():
			line = []
			slope_info = cmpd.curve_data.get_slope_CI(CI_val = self.options['LC_CI'])
			rel_pot = self.compare_LC(cmpd1 = cmpd_name, 
				cmpd2 = self.options['REFERENCE_COMPOUND'], 
				LCval = self.options['LC_VALUES'],
				CI = self.options['REL_POT_CI'],
				n_bs = 10000)

			rel_done = False
			LC_done = False

			for item in header:
				if item.lower() in 'compound': line.append(cmpd.data["name"])
				elif 'trial' in item.lower(): line.append(str(cmpd.data["n_trials"]))
				elif item.lower() in 'slope': line.append(format_lc_val(slope_info[1],3))
				elif 'slope' in item.lower() and 'ci' in item.lower(): 
					line.append(CI_to_string(slope_info[0], slope_info[2]))
				elif self.options['REFERENCE_COMPOUND'].lower() in item.lower():
					if rel_done: continue
					else:
						line = [*line, *format_LC_to_CSV(rel_pot)]
						rel_done = True
				elif self.options['REFERENCE_COMPOUND'].lower() not in item.lower() and "lc" in item.lower():
					if LC_done: continue
					else:
						line = [*line, *format_LC_to_CSV(cmpd.get_LC_CIs(**self.options))]
						LC_done = True
				elif "codes" in item.lower(): line.append('"' + ",".join([i for i in cmpd.data["ids"]]) + '"')
				elif "date" in item.lower(): line.append('"' + ",".join([i for i in cmpd.data["test_dates"]]) + '"')
				else: line.append(" ")
			output.append(",".join(line))
		return output

	def generate_csv_header(self):
		LC_title_names = []
		for idx, LC in enumerate(self.options['LC_VALUES']): 
			LC_title_names.append("LC" + str(round(self.options['LC_VALUES'][idx] * 100)))
			LC_title_names.append("LC" + str(round(self.options['LC_VALUES'][idx]* 100)) + 
						" " + str(round(self.options['LC_CI']* 100)) + "%CI")
		rel_pot_col = []
		for idx, LC in enumerate(self.options['LC_VALUES']): 
			rel_pot_col.append(self.options['REFERENCE_COMPOUND'] + ' Pot. Rel. to Cmpd. at LC' + 
						str(round(self.options['LC_VALUES'][idx] * 100)))
			rel_pot_col.append(self.options['REFERENCE_COMPOUND'] + ' Pot. Rel. to Cmpd. at LC' + 
						str(round(self.options['LC_VALUES'][idx] * 100)) + 
						" " + str(round(self.options['REL_POT_CI']* 100)) + "%CI")
		header = ['Compound', 
					'Number of Trials', 
					*LC_title_names,
		 			'Slope', 
		 			'Slope CI',  
		 			*rel_pot_col,
		 			'Tested Codes', 
		 			'Test Dates']
		return header

	def save_plots(self, name, image_dir, image, *args, **kwargs):
		return

	def compare_LC(self, cmpd1, cmpd2, LCval = np.array([0.5]), CI = 0.95, n_bs = 10000, *args, **kwargs):
		'''
		Calculates relative potency using random samples from the kernel density of the LCx distribution for
		cmpd1 and cmpd2.
		'''
		kernel1 = self.cmpd_data[cmpd1].curve_data.LC_kernel(LCval).sample(n_samples=n_bs)
		kernel2 = self.cmpd_data[cmpd2].curve_data.LC_kernel(LCval).sample(n_samples=n_bs)
		diff = kernel1 - kernel2
		lb = (1. - CI)/2.
		ub = 1.-lb
		quantiles = [lb, 0.5, ub]
		vals = np.quantile(diff, quantiles, interpolation='linear', axis = 0)
		return np.exp(vals).T

	def full_process(self, 
						new_datafile = None, 
						key_file = None,
						csv_outfile = None, 
						archive_path = None,
						image_dir = None, 
						*args, 
						**kwargs):
		self.merge_old_new(new_datafile, archive_path, key_file)
		# self.process_compounds(*args, **kwargs)

		output_info = []

		for cmpd_name, cmpd in self.cmpd_data.items():
			tic = time()
			print(cmpd.data["name"])
			cmpd.fit_data(options = self.options)
			ax = cmpd.make_plot()
			if image_dir: self.save_plots(name, image_dir, ax)
 
		header = self.generate_csv_header()
		body = self.generate_csv_data_lines(header)
			
		if csv_outfile: self.save_csv(csv_outfile, header, body, *args, **kwargs)
		if archive_path: self.save_archive(archive_path, *args, **kwargs)
		

if __name__ == "__main__":
	MA = MerlinAnalyzer()
	MA.full_process(new_datafile = "./test.csv", 
		key_file = "./key.csv", 
		archive_path=".", 
		csv_outfile = './output.csv', 
		image_dir= "./images")


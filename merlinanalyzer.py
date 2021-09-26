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
from latex_writer import LatexWriter

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
		cmpd_data.dropna(subset = ['Live', 'Dead', 'Compound'], inplace=True)
		cmpd_data = cmpd_data[cmpd_data.Count != 0]
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
			for cmpd in self.cmpd_data.keys():
				dict_change = utils.check_library_change(self.cmpd_data[cmpd].options, self.options)
				if dict_change: self.cmpd_data[cmpd].curve_data.reset_curves()
				#update dictionaries
				for k, v in self.options.items():
					self.cmpd_data[cmpd].options[k] = v
					self.cmpd_data[cmpd].curve_data.options[k] = v


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
							min_conc = min(cmpd_data["Conc"].tolist()),
							n_trials = len(set(unique_ids)),
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
		output = []

		for cmpd_name, cmpd in self.cmpd_data.items():
			line = []
			comment = " "
			good_curve = True
			#Check to make sure that the slope is positive enough
			slope_info = cmpd.curve_data.get_slope_CI(CI_val = 0.95)

			if slope_info[1] < 1.e-2:
				comment += "Fitted slope is too shallow. "
				good_curve = False

			#Check to make sure that the LC50 value is close enough to the 
			LC50_info = np.power(2.,cmpd.curve_data.get_LC50_CI(CI_val=0.95, log = True))
			lc_vals = utils.format_LC_to_CSV(cmpd.get_LC_CIs()) 
			# print(cmpd_name, LC50_info[1], 2**(1 + cmpd.data["max_conc"]) ,  2**(cmpd.data["min_conc"] - 1))

			if LC50_info[1] > 2**(1 + cmpd.data["max_conc"]) or LC50_info[1]< 2**(cmpd.data["min_conc"] - 1) : 
				# print(LC50_info[1], {lc_vals[0]}, cmpd.data["min_conc"]/2., 2*cmpd.data["max_conc"])
				comment += f"Calculated LC50 ({utils.format_lc_val(LC50_info[1])}) out of bounds. "
				lc_vals = ['NA'] * len(lc_vals)
				good_curve = False
				
			if good_curve:
				rel_pot = self.compare_LC(cmpd = cmpd_name, n_bs = 100000)
				rel_pot = utils.format_LC_to_CSV(rel_pot)
			else: 
				lc_vals = ['NA'] * len(lc_vals)
				rel_pot = ["NA"] * 2*len(self.options['LC_VALUES'])

			rel_done = False
			LC_done = False

			for item in header:
				if item.lower() in 'compound': line.append(cmpd.data["name"])
				elif 'rows' in item.lower(): line.append(str(cmpd.data["n_trials"]))
				elif item.lower() in 'slope': 
					line.append(utils.format_lc_val(slope_info[1]))
				elif 'slope' in item.lower() and 'ci' in item.lower(): 
					line.append(utils.CI_to_string(slope_info[0], slope_info[2]))
				elif self.options['REFERENCE_COMPOUND'].lower() in item.lower():
					if rel_done: continue
					else:
						line = [*line, *rel_pot]
						rel_done = True
				elif self.options['REFERENCE_COMPOUND'].lower() not in item.lower() and "lc" in item.lower():
					if LC_done: continue
					else:
						line = [*line, *lc_vals]
						LC_done = True
				elif "codes" in item.lower(): line.append('"' + ", ".join(list(set([i for i in cmpd.data["ids"]]))) + '"')
				elif "date" in item.lower(): line.append('"' + ", ".join(list(set([i for i in cmpd.data["test_dates"]]))) + '"')
				elif "bio" in item.lower(): line.append(f"{len(cmpd.data['test_dates'])}")
				elif item == "R2": line.append(utils.format_lc_val(cmpd.curve_data.r2))
				elif "comment" in item.lower():
					comment = comment if len(comment) == 1 else comment[1:len(comment)] 
					line.append(comment)
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
		if self.options['REL_POT_TO_REF']:
			init_name = self.options['REFERENCE_COMPOUND'] + ' Pot. Rel. to Cmpd at LC'
		else:
			init_name = "Pot. of Cmpd Rel. to " + self.options['REFERENCE_COMPOUND'] + " at LC"
		for idx, LC in enumerate(self.options['LC_VALUES']): 
			rel_pot_col.append(init_name + str(round(self.options['LC_VALUES'][idx] * 100)))
			rel_pot_col.append(init_name + str(round(self.options['LC_VALUES'][idx] * 100)) + 
						" " + str(round(self.options['REL_POT_CI']* 100)) + "%CI")
		header = ['Compound', 
					'Biological Reps',
					'Total Rows', 
					*LC_title_names,
					*rel_pot_col,
		 			'R2',
		 			'Slope', 
		 			'Slope CI',  
		 			'Tested Codes', 
		 			'Test Dates',
		 			'Comments']
		return header

	def compare_LC(self, cmpd, n_bs = 100000, *args, **kwargs):
		'''
		Calculates relative potency using random samples from the kernel density of the LCx distribution for
		cmpd1 and cmpd2.
		'''
		cmpd_kernel = self.cmpd_data[cmpd].curve_data.LC_kernel(
			LC_val = 1. - self.options['LC_VALUES']).sample(n_samples=n_bs)
		ref_kernel= self.cmpd_data[self.options['REFERENCE_COMPOUND']].curve_data.LC_kernel(
			LC_val = 1. - self.options['LC_VALUES']).sample(n_samples=n_bs)
		diff = cmpd_kernel - ref_kernel if self.options['REL_POT_TO_REF'] else kernel2 - kernel1
		func = utils.calc_ET_CI if self.options["CI_METHOD"].lower() in utils.ET_VARS else utils.calc_HPDI_CI
		vals = func(diff, CI_level = self.options['REL_POT_CI'])

		return np.power(2., vals).T

	def full_process(self, 
						new_datafile = None, 
						key_file = None,
						csv_outfile = None, 
						out_path = None,
						pdf_outfile = 'graphs.pdf',
						archive_path = None,
						*args, 
						**kwargs):
		self.merge_old_new(new_datafile, archive_path, key_file)
		# self.process_compounds(*args, **kwargs)

		image_dir = os.path.abspath(os.path.join(out_path, 'images'))
		pdf_dir = os.path.abspath(os.path.join(image_dir, 'pdf'))

		if not os.path.exists(image_dir):
			os.makedirs(image_dir)
		if not os.path.exists(pdf_dir):
			os.makedirs(pdf_dir)

		LW = LatexWriter(img_folder = pdf_dir)
		# cmpd_ct = 0
		for cmpd_name, cmpd in self.cmpd_data.items():
			cmpd.fit_data(options = self.options)
			cmpd.make_plot()
			pdf_path = cmpd.plot.save_plot(cmpd_name, image_dir, pdf_dir = pdf_dir)
			# print(cmpd.curve_data.get_LC50_CI(CI_val=0.95).squeeze())
			lc50lb, lc50med, lc50ub = cmpd.curve_data.get_LC50_CI(CI_val=0.95).squeeze()

			if lc50med > 2**(1 + cmpd.data["max_conc"]) or lc50med< 2**(cmpd.data["min_conc"] - 1): lcTrue = False 
			else: lcTrue = True 
				
			lc50med = utils.format_lc_val(lc50med)
			lc50CI = '[' + utils.format_lc_val(lc50lb) + ', ' + utils.format_lc_val(lc50ub) + ']'
			LW.make_cmpd_graph(image_dir = pdf_path, 
				name = cmpd_name, 
				lc50 = lc50med, 
				lc50CI = lc50CI, 
				lcTrue = lcTrue,
				R2 = utils.format_lc_val(cmpd.curve_data.r2), 
				reps = len(cmpd.data['test_dates']))
		self.save_archive(archive_path, *args, **kwargs)
		header = self.generate_csv_header()
		body = self.generate_csv_data_lines(header)
		if csv_outfile: self.save_csv(os.path.abspath(os.path.join(out_path, csv_outfile)), header, body, *args, **kwargs)
		LW.write_file(out_path = os.path.abspath(os.path.join(out_path, pdf_outfile+".tex")) )
		# LW.make()

import numpy as np
from curvell import CI_finder
import matplotlib.pyplot as plt
from compound import Compound

class MerlinAnalyzer:

	def __init__(self, filename, *args, **kwargs):
		self.filename = filename
		self.data = {}

	def worm_datafile_check(self, *args, *kwargs):
		def inner(*args, **kwargs):
			try:
				file = open(self.filename, 'r')
				close(file)
			except FileNotFoundError: 
			#TODO: check for proper worm formatting
			func()
			return
		return

	@worm_datafile_check
	def read_new_data(self, filename):
		return

	@worm_datafile_check
	def read_old_data(self, json_filename):
		return

	def save_json(self, json_filename):
		return

	def save_csv(Self, csv_filename):
		return

	def compare_LC(self, cmpd1, cmpd2, LCval = 50, CI = 0.95):
		return
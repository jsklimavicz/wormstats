import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import os

class MerlinGrapher:
	'''
	Class for holding and modifying graphs for data and dost-response data.
	'''
	def __init__(self, *args, options = {}, **kwargs):
		self.f, self.ax = plt.subplots()
		self.data = kwargs
		self.options = options

	def plot_data(self):
		'''
		Plots the actual data.
		'''
		#Confidence intervals for the curve.
		if "x" in self.data and "lb" in self.data and "ub" in self.data and self.options["PLOT_CURVE_ERROR"]:
			self.ax.fill_between(self.data["x"], 
				self.data["lb"], 
				self.data["ub"], 
				color = self.options["CURVE_CI_COLOR"],
				alpha = self.options["ALPHA"] if "ALPHA" in self.options else 0.5)
		#Data points.
		if "conc" in self.data and "probs" in self.data and self.options["PLOT_DATA_POINTS"]:
			conc = self.data["conc"].copy()
			if self.options["JITTER"]: conc += np.random.uniform(-self.options["JITTER_FACTOR"], self.options["JITTER_FACTOR"], len(self.data["conc"]))
			self.ax.plot(conc, 
				self.data["probs"], 
				marker = self.options["MARKER_STYLE"], 
				mew = 0.0, 
				mfc = self.options["POINT_COLOR"], 
				ls = 'None')
			#Error bars. 
			if "error_bars" in self.data and self.options["PLOT_ERROR_BARS"] and self.options["ERROR_BAR_CUTOFF"] > self.data["n_trials"]:
				self.ax.errorbar(x = conc, 
					y = self.data["probs"], 
					yerr = self.data["error_bars"], 
					ecolor = self.options["POINT_COLOR"],
					elinewidth=0.5,
					ls = 'None')
		#Best fit line.
		if "x" in self.data and "line" in self.data and self.options["PLOT_LINE"]: 
			self.ax.plot(self.data["x"], 
				self.data["line"], 
				c = self.options["LINE_COLOR"], 
				ls = self.options["LINE_STYLE"])

	def set_labels(self, title):
		'''
		Set axes labels, ticks, title, etc. 
		'''
		self.ax.set_xlabel('Concentration (ppm)')
		self.ax.set_title(label = title)
		self.ax.set_ylabel('Percent Survival')
		xticks, xticklabels = self.calc_x_ticks()
		plt.xticks(ticks=xticks, rotation=90, fontsize=8, labels=xticklabels)
		plt.yticks(ticks=np.array(range(0, 5, 1))/4. , labels=np.array(range(0, 101, 25)) )
		plt.grid(b=True, alpha = 0.5)
		plt.xlim([min(xticks), max(xticks)])
		plt.ylim([0,1])
		plt.gcf().subplots_adjust(bottom=0.15)

	def calc_x_ticks(self):
		'''
		Calculate the x-ticks for the graph based on the range of concentrations in the data. 
		'''
		lb, ub = round(min(self.data["conc"])), round(max(self.data["conc"]))
		xticks = np.array(range(lb-1, ub+2, 1))
		xticklabels = [round(2**i) if i >=0 else 2.**i for i in xticks]
		for i, label in enumerate(xticklabels):
			if label < 1./64.: xticklabels[i] = "{:.1e}".format(label).replace("e-0", "e-")
			elif label < 1./4.: xticklabels[i] = "{:.3f}".format(label)
			else: xticklabels[i] = str(label) 
		return xticks, xticklabels

		
	def save_plot(self, name, image_dir, *args, pdf = True, pdf_dir = None, close_after = True, **kwargs):
		filename = name + ".png"
		png_path = os.path.join(image_dir, filename)
		plt.savefig(png_path)
		if pdf:
			if pdf_dir is None: pdf_dir = os.path.join(image_dir, "pdf")
			pdf_path = os.path.join(pdf_dir, name + ".pdf")
			plt.savefig(pdf_path, format='pdf',bbox_inches='tight', transparent=True)
		if close_after: plt.close()
		return pdf_path

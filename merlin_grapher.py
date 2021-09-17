import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

class MerlinGrapher:
	def __init__(self, *args, options = {}, **kwargs):
		self.f, self.ax = plt.subplots()
		self.data = kwargs
		self.options = options

	def plot_data(self):
		if "x" in self.data and "lb" in self.data and "ub" in self.data:
			self.ax.fill_between(self.data["x"], 
				self.data["lb"], 
				self.data["ub"], 
				color = self.options["CURVE_CI_COLOR"],
				alpha = self.options["ALPHA"] if "ALPHA" in self.options else 0.5)
		if "x" in self.data and "line" in self.data: 
			self.ax.plot(self.data["x"], 
				self.data["line"], 
				c = self.options["LINE_COLOR"], 
				ls = self.options["LINE_STYLE"])
		if "conc" in self.data and "probs" in self.data:
			conc = self.data["conc"]
			print(conc)
			if self.options["JITTER"]: conc += np.random.uniform(-self.options["JITTER_FACTOR"], self.options["JITTER_FACTOR"], len(self.data["conc"]))
			print(conc)
			self.ax.plot(conc, 
				self.data["probs"], 
				marker = self.options["MARKER_STYLE"], 
				mew = 0.0, 
				mfc = self.options["POINT_COLOR"], 
				ls = 'None')
			if "error_bars" in self.data:
				self.ax.errorbar(x = conc, 
					y = self.data["probs"], 
					yerr = self.data["error_bars"], 
					ecolor = self.options["POINT_COLOR"],
					ls = 'None')


	def set_labels(self, title):
		self.ax.set_xlabel('Concentration (ppm)')
		self.ax.set_title(label = title)
		self.ax.set_ylabel('Percent Survival')
		xticks, xticklabels = self.calc_x_ticks()
		plt.xticks(ticks=xticks, rotation=90, fontsize=8, labels=xticklabels)
		plt.yticks(ticks=np.array(range(0, 5, 1))/4. , labels=np.array(range(0, 101, 25)) )
		plt.grid(b=True, alpha = 0.5)
		plt.xlim([min(xticks), max(xticks)])
		plt.ylim([0,1])

	def calc_x_ticks(self):
		lb, ub = round(min(self.data["conc"])), round(max(self.data["conc"]))
		xticks = np.array(range(lb-1, ub+2, 1))
		xticklabels = [round(2**i) if i >=0 else 2.**i for i in xticks]
		for i, label in enumerate(xticklabels):
			if label < 0.125: xticklabels[i] = "{:.2e}".format(label)
			else: xticklabels[i] = str(label) 
		return xticks, xticklabels
		



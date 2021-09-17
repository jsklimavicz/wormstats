import matplotlib.pyplot as plt
import numpy as np

class MerlinGrapher:
	def __init__(self, *args, **kwargs):
		self.f, self.ax = plt.subplots()
		self.data = kwargs

	def plot_data(self, *args, **kwargs):
		for k, v in kwargs.items(): self.data[k] = v
		if "x" in self.data and "lb" in self.data and "ub" in self.data:
			self.ax.fill_between(self.data["x"], self.data["lb"], self.data["ub"], alpha = self.data["alpha"] if "alpha" in self.data else 0.5)
		if "x" in self.data and "line" in self.data: 
			self.ax.plot(self.data["x"], self.data["line"], 'r-')
		if "conc" in self.data and "probs" in self.data:
			self.ax.plot(self.data["conc"], self.data["probs"], 'g.')

	def set_labels(self, title):
		self.ax.set_xlabel('Concentration (ppm)')
		self.ax.set_title(label = title)
		self.ax.set_ylabel('Percent Survival')
		xticks, xticklabels = self.calc_x_ticks()
		self.ax.set_xticks(ticks=xticks)
		self.ax.set_xticklabels(labels=xticklabels)
		self.ax.set_yticks(ticks=np.array(range(0, 5, 1))/4. )
		self.ax.set_yticklabels(labels=np.array(range(0, 101, 25)) )
		plt.grid(b=True, alpha = 0.5)
		plt.xlim([min(xticks), max(xticks)])
		plt.ylim([0,1])

	def calc_x_ticks(self):
		lb, ub = round(min(self.data["conc"])), round(max(self.data["conc"]))
		xticks = np.array(range(lb-1, ub+1, 1))
		xticklabels = [round(2**i) if i >=0 else 2.**i for i in xticks]
		return xticks, xticklabels
		



#!/usr/bin/env python3
#late_writer.py
import sys
import argparse
import os, errno
import glob
import re
import time
import numpy as np
import math
import pandas as pd
import time



class Latex_Writer:
	gridx = 2
	gridy = 4
	bs = '\\'

	def __init__(self, title = "Merlin Bioassay Results"):
		self.page_max = args.p
		self.make_header()
		self.make_title(title)
		self.cmpd_graph_list = []
		self.graph_count = 0
		# self.section_driver()
		# self.make_end()
		# self.write_file()


	def make_header(self):
		self.header = ["\\documentclass{article}"]
		# self.header.append("\\usepackage[utf8]{inputenc}")
		self.header.append("\\usepackage[letterpaper,left=0.75in, right=0.75in, top=0.75in, bottom=0.75in]{geometry}")
		#Packages for imagesPackages for images
		self.header.append("\\usepackage{graphicx}")
		# self.header.append("\\usepackage[dvipsnames]{xcolor}")
		# self.header.append("\\definecolor{TTgreen}{cmyk}{0.56, 0, 82, 0.26}")
		# self.header.append("\\usepackage{svg}")
		self.header.append("\\usepackage[font=scriptsize]{subcaption}")
		self.header.append("\\captionsetup[subfigure]{format=hang,justification=raggedright,singlelinecheck=false}")
		#Packages for symbols and chemistry text
		self.header.append("\\usepackage{amsmath,amsthm,amsfonts,amssymb,mathtools}")
		# self.header.append("\\usepackage[version=4]{mhchem}")
		self.header.append(f"\\graphicspath{{{{{self.image_folder:s}}}}}")
		self.header.append("\\setcounter{topnumber}{8}")
		self.header.append("\\setcounter{bottomnumber}{8}")
		self.header.append("\\setcounter{totalnumber}{8}")
		self.header.append("\\usepackage{fancyhdr}")
		self.header.append("\\pagestyle{fancy}")
		self.header.append("\\fancyhf{}")
		self.header.append("\\fancyhead[L]{\\rightmark}")
		self.header.append("\\fancyfoot[C]{\\thepage}")
		self.header.append("\\fancyhead[R]{\\today}")
		self.header.append("\\renewcommand{\\headrulewidth}{0pt}")
		self.header.append("\n\n")


		self.header.append("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")


	def make_title(self, title):
		self.title = [f"\\title{{{title:s}}}"]
		self.title.append("\\author{James Klimavicz}")
		self.title.append("\\date{\\today}")
		self.title.append("\\begin{document}")
		self.title.append("\n\n")


	def make_end(self):
		self.end = ["\n\n"]
		self.end.append("\\end{document}")

	def make_cmpd_caption(self, name, lc50, lc50CI, R2, reps):
		'''
		<Compound Name>
		LC50: <LC50> <CI>; Curve r^2: <R2> 
		Biological Replicates: <Bio reps>
		'''
		caption = f'{name} \\\\'
		caption += 'LC_{50} ' + f'{lc50} ppm {lc50CI} \\\\'
		if reps == 1:
			caption += f'{reps} biologicial replicate; R^2: {r2}'
		else:
			caption += f'{reps} biologicial replicates; R^2: {r2}'
		#TODO: add reference compound rel. potency
		return caption

	def make_cmpd_graph(self, img_path, name, lc50, lc50CI, R2, reps):
		graph_lines = []
		img_line = f"\\includegraphics[width = {{0.95\\textwidth}}]{{{img_path:s}}}"
		name, alpha_name = latexify_name(name)
		caption = make_cmpd_caption(name, lc50, lc50CI, R2, reps)
		graph_lines.append(f"   \\begin{{subfigure}}{{{image_width:.3f}\\textwidth}}")
		graph_lines.append("      \\centering")
		graph_lines.append("      "+img_line)
		graph_lines.append("      \\vspace{-0.35cm}")
		graph_lines.append("      "+caption)
		graph_lines.append("   \\end{subfigure}%")
		self.graph_count += 1
		self.cmpd_graph_list.append(alpha_name, graph_lines) 


	def make_body(self):
		self.body = []
		img_count = 0
		pages = 0
		row_end = False
		image_width = 1/self.gridx
		self.cmpd_graph_list = self.cmpd_graph_list.sort(key=lambda x:x[0])
		for index, (alpha_name, graph_lines) in enumerate(self.cmpd_graph_list.items()):
			if img_count % (self.gridx * self.gridy) == 0:  self.body.append("\\begin{figure}[thp!]")
			self.body.append(i) for i in graph_lines
			self.body.append("\n")
			img_count += 1
			if img_count % (self.gridx * self.gridy) == 0: 
				self.body.append("\\end{figure}")
				img_count = 0
				pages += 1
				if pages >= self.page_max:
					break
				self.body.append("\\clearpage")
				self.body.append("\\pagebreak")
			if img_count % self.gridx == 0: self.body.append("\\vspace{-0.1cm}")

		if img_count % (self.gridx * self.gridy) != 0:
			self.body.append("\\end{figure}")
		self.body.append("\\pagebreak")
		self.body.append("\n\n")



	def write_file(self):
		file = open(self.output, 'w')
		file.write("\n".join(self.header))
		file.write("\n".join(self.title))
		file.write("\n".join(self.body))
		file.write("\n".join(self.end))
		file.close()


def latexify_name(name):
	greek_letters = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'lambda']
	mod = False
	alpha_name = name
	if name[:2].lower() == "s-": 
		alpha_name = name[2:]
		name = "\\textit{S}-" + alpha_name
		mod = True
	elif name[:2].lower() == "r-": 
		alpha_name = name[2:]
		name = "\\textit{R}-" + alpha_name
		mod = True
	else:
		for letter in greek_letters:
			if name[:len(letter)+1].lower() = letter+'-': 
				alpha_name = name[len(letter):]
				name = "\\" + letter + alpha_name
				mod = True
				break
	if len(name)< 5:
		return name.upper(), alpha_name.lower()
	elif mod:
		return name, alpha_name.lower()
	else:
		return name[0].capitalize() + name[1:], alpha_name.lower()

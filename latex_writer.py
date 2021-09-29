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



class LatexWriter:
	'''
	Class for making a LaTeX file and corresponding .pdf from all of the graphed data. 
	'''

	#control number of plots on a page.
	gridx = 2
	gridy = 3

	def __init__(self, img_folder, title = "Merlin Bioassay Results"):
		self.image_folder = img_folder
		self.make_header()
		self.make_title(title)
		self.cmpd_graph_list = []
		self.graph_count = 0
		self.make_end()

	def make_header(self):
		'''
		Header section for the LaTeX file.
		'''
		self.header = ["\\documentclass{article}"]
		self.header.append("\\usepackage[letterpaper,left=0.75in, right=0.75in, top=0.75in, bottom=0.75in]{geometry}")
		self.header.append("\\usepackage{graphicx}")
		self.header.append("\\usepackage[justification=centering]{subcaption}")
		self.header.append("\\captionsetup[subfigure]{format=hang,justification=raggedright,singlelinecheck=false}")
		self.header.append("\\usepackage{amsmath,amsthm,amsfonts,amssymb,mathtools}")
		self.header.append("\\usepackage[mmddyyyy,HHmmss]{datetime}")
		self.header.append(f"\\graphicspath{{{{{self.image_folder:s}}}}}")
		self.header.append("\\setcounter{topnumber}{8}")
		self.header.append("\\setcounter{bottomnumber}{8}")
		self.header.append("\\setcounter{totalnumber}{8}")
		self.header.append("\\usepackage{fancyhdr}")
		self.header.append("\\pagestyle{fancy}")
		self.header.append("\\fancyhf{}")
		self.header.append("\\fancyhead[L]{Merlin Bioassay Results}")
		self.header.append("\\fancyfoot[C]{\\thepage}")
		self.header.append("\\fancyhead[R]{Compiled on \\today\\ at \\currenttime}")
		self.header.append("\\renewcommand{\\headrulewidth}{0pt}")
		self.header.append("\n\n")
		self.header.append("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")


	def make_title(self, title):
		'''
		Title section of document 
		'''
		self.title = [f"\\title{{{title:s}}}"]
		self.title.append("\\author{James Klimavicz}")
		self.title.append("\\date{\\today}")
		self.title.append("\\begin{document}")
		self.title.append("\n\n")


	def make_end(self):
		'''
		Portion to close out the LaTeX document. 
		'''
		self.end = ["\n\n"]
		self.end.append("\\end{document}")

	def make_cmpd_caption(self, name, lc50, lc50CI, R2, reps, lcTrue = True):
		'''
		Creates the caption for each image. Retruns the caption as a string. 
		**<Compound Name>** LC50: <LC50> ppm <CI>
		<Bio reps> biological replicates; R^2: <R2> 
		'''
		caption = f'\\textbf{{{name}}}' 
		if lcTrue: caption += ' LC$_{50}$: ' + f'{lc50} ppm {lc50CI} \\\\ \n'
		else: caption += ' LC$_{50}$ cannot be estimated with the given data. \\\\ \n'

		if reps == 1: caption += f'{reps} biological replicate; R$^2$: {R2}'
		else: caption += f'{reps} biological replicates; R$^2$: {R2}'
		#TODO: add reference compound rel. potency
		return caption

	def make_cmpd_graph(self, image_dir, name, lc50, lc50CI, R2, reps, lcTrue = True):
		'''
		Portion that the graphic and caption for a given compound. Input values are 
		from the compound info from curve fitting. 
		'''
		graph_lines = []
		image_width = 1/self.gridx
		img_line = f"\\includegraphics[width = {{0.95\\textwidth}}]{{{image_dir:s}}}"
		name, alpha_name = self.latexify_name(name)
		print(alpha_name, name)
		caption = f'\\caption*{{{self.make_cmpd_caption(name, lc50, lc50CI, R2, reps, lcTrue = lcTrue)}}}'
		graph_lines.append(f"   \\begin{{subfigure}}{{{image_width:.3f}\\textwidth}}")
		graph_lines.append("      \\centering")
		graph_lines.append("      "+img_line)
		graph_lines.append("      \\vspace{-0.05cm}")
		graph_lines.append("      "+caption)
		graph_lines.append("      \\vspace{0.1cm}")
		graph_lines.append("   \\end{subfigure}%")
		self.graph_count += 1
		self.cmpd_graph_list.append((alpha_name, graph_lines))


	def make_body(self):
		'''
		Driver for including all figures/captions for the output file. 
		'''
		self.body = [] #array of lines. 
		self.graph_count = 0
		row_end = False

		#sort by alphabetizing name, determined by latexify_name
		self.cmpd_graph_list.sort(key=lambda x:x[0]) 

		#include a figure of subfigures with subcaptions. Layout determined by self.gridx and self.gridy
		for index, (alpha_name, graph_lines) in enumerate(self.cmpd_graph_list):
			if self.graph_count % (self.gridx * self.gridy) == 0:  self.body.append("\\begin{figure}[thp!]")
			for i in graph_lines: self.body.append(i)
			self.graph_count += 1
			if self.graph_count % (self.gridx * self.gridy) == 0: 
				self.body.append("\\end{figure}")
				self.graph_count = 0
				self.body.append("\\clearpage")
				self.body.append("\\pagebreak")
			if self.graph_count % self.gridx == 0: self.body.append("\\vspace{-0.1cm}")
		if self.graph_count % (self.gridx * self.gridy) != 0: self.body.append("\\end{figure}")
		self.body.append("\\pagebreak")
		self.body.append("\n\n")

	def write_file(self, out_path):
		'''
		Creates the actual LaTeX file at out_path
		'''
		self.make_body()
		file = open(out_path, 'w')
		file.write("\n".join(self.header))
		file.write("\n".join(self.title))
		file.write("\n".join(self.body))
		file.write("\n".join(self.end))
		file.close()
		self.make(out_path)

	def make(self, out_path):
		'''
		Convert the LaTeX file out_path to a pdf in the same directory. 
		'''
		pwd = os.getcwd()
		os.chdir(os.path.dirname(out_path))
		os.system(f'pdflatex -interaction=nonstopmode {out_path:s} --shell-escape  --enable-write18  > pdf_tex.log 2>&1')
		os.chdir(pwd)

	def latexify_name(self, name):
		'''
		Stylizes the name of a compound to include LaTeX-style letters, e.g. greek letters, italicized stereochem, etc.
		Currently, R- and S- stereochem designators and lowercase greek letters are implemented. Only looks at the 
		start of the name to ensure that middles of names are not replaced, though this behavior may not be desirable later. 
		Returns the LaTeX-style name and a name without these designators for the purpose of alphabetizing.
		'''
		greek_letters = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta', 'iota',
							'kappa', 'lambda', 'mu', 'nu','xi', 'omicron','pi', 'rho', 'sigma','tau', 'upsilon',
							'phi', 'chi','psi','omega']
		mod = False
		alpha_name = name

		#stereochem designators
		if name[:2].lower() == "s-": 
			alpha_name = name[2:]
			name = "\\textit{S}-" + alpha_name
			mod = True
		elif name[:2].lower() == "r-": 
			alpha_name = name[2:]
			name = "\\textit{R}-" + alpha_name
			mod = True
		else:
			#Greek letters
			for letter in greek_letters:
				if name[:len(letter)+1].lower() == letter+'-': 
					alpha_name = name[len(letter)+1:]
					name = "$\\" + letter + '$-' + alpha_name.capitalize()
					mod = True
					break
		#Return name and alphabetizing name 
		if len(name)< 5:
			return name.upper(), alpha_name.lower()
		elif mod:
			return name, alpha_name.lower()
		else:
			return name[0].capitalize() + name[1:], alpha_name.lower()

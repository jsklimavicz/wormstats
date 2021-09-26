from merlinanalyzer import MerlinAnalyzer
import sys
import argparse
import os
import errno

parser = argparse.ArgumentParser(description="""Produces pdf of figures.""")
parser.add_argument("-i", "--input", metavar="/path/to/input.csv", required = True,
	help="Input live/count file for Merlin in .csv format.")
parser.add_argument("-k", "--key", metavar="/path/to/key.csv", required = True,
	help="Key file that contains the compound codes and names.")
parser.add_argument("-p", "--path", metavar="/path/for/saving", required = True,
	help="Directory in which the output will be saved.")
parser.add_argument("-o", "--output_csv", metavar="output_filename.csv", default = "output.csv",
	help="Name of the output .csv file.")
parser.add_argument("-a", "--archive_path", metavar="/path/to/pickle/archive", default = '.',
	help="Name of the output .csv file.")
parser.add_argument("-g", "--output_pdf", metavar="graphs_filename", default = "graphs",
	help="Name of the output file containing the graphs WITHOUT file extension.")
# parser.add_argument("-o", "--output", metavar="output_filename.csv", default = "output.csv",
# 	help="Name of the output .csv file.")


if __name__ == "__main__":
	args = parser.parse_args()
	args.archive_path = args.path
	MA = MerlinAnalyzer()
	MA.full_process(new_datafile = args.input, 
		key_file = args.key, 
		out_path=args.path, 
		archive_path = args.archive_path,
		csv_outfile = args.output_csv,
		pdf_outfile = args.output_pdf)

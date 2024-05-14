#!/usr/bin/env python3
from os import listdir
from os.path import isfile, join
import sys
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy.interpolate import interp1d


if __name__ == "__main__":

	path = "./"
	save = 1
	file_print = None
	files = [f for f in listdir(path) if isfile(join(path, f))]
	energy = sys.argv[1] 
	for file in files:
		if (float(energy) > 0):
			if ("positions_" in file) and (energy in file) and ("all" not in file) and ("-" not in file):
			# if ("positions_all" in file) and (energy in file) and ("-" not in file):
				file_print = file
				break
		else:
			if ("positions_" in file) and (energy in file) and ("all" not in file):
			# if ("positions_all" in file) and (energy in file):		
				file_print = file
				break

	if (file_print == None):
		raise ValueError("No file for selected energy")

	with open(file_print, "r") as file:
		lines = file.readlines()

	size = int(( len(lines)**0.5) / 0.4 )	# side of "particle" square equal to 0.4 of total side
	print("Droplet's size is {}".format(size))

	if 1:
		arr = np.zeros((size, size))
		x_coords = []
		y_coords = [] 
		for line in lines:
			line_c = line.split()
			x = int(line_c[0])
			y = int(line_c[1])
			arr[x][y] = 1

		if 1:
			fig, ax = plt.subplots()
			ax.tick_params(axis='both', which='major', labelsize=16)
			ax.tick_params(axis='both', which='minor', labelsize=12)
			plt.imshow(arr)
			# plt.show()
			file_save = file_print.replace("txt", "png")
			if (save):
				plt.savefig("figures/"+file_save,bbox_inches='tight', dpi=600)

	
	file_print = None
	files = [f for f in listdir(path) if isfile(join(path, f))]
	energy = sys.argv[1] 
	for file in files:
		# if ("energies_all" in file) and (energy in file):
		if (float(energy) > 0):
			if ("energies_" in file) and (energy in file) and ("all" not in file) and ("-" not in file):
				file_print = file
				break
		else:
			if ("energies_" in file) and (energy in file) and ("all" not in file):
				file_print = file
				break

	if (file_print == None):
		raise ValueError("No file for selected energy")

	with open(file_print, "r") as file:
		lines = file.readlines()

	eff_list = []
	efs_list = []
	step_list = []
	for line in lines:
		line_c = line.split()
		eff_list.append(float(line_c[0]))
		efs_list.append(float(line_c[1]))
		step_list.append(float(line_c[2]))


	# Efs plot 
	fig, ax = plt.subplots()
	fig.set_size_inches(9, 9)

	ax.set_ylabel(r'Energy per particle (kT)', fontsize=30)
	ax.set_xlabel(r'Time (step)', fontsize=30)
	ax.tick_params(axis='both', which='major', labelsize=26)
	ax.tick_params(axis='both', which='minor', labelsize=20)
	ax.xaxis.offsetText.set_fontsize(22)
	plt.gcf().subplots_adjust(left=0.15)

	# ax.plot(step_list, eff_list, '-', c='red', label="Eff")
	lowess = sm.nonparametric.lowess(efs_list, step_list, frac=.09)
	lowess_x = list(zip(*lowess))[0]
	lowess_y = list(zip(*lowess))[1]
	f = interp1d(lowess_x, lowess_y, bounds_error=False)
	xnew = np.linspace(step_list[0], step_list[-1], 10000)
	ynew = f(xnew)
	ax.plot(xnew, ynew, '-', c='red', linewidth=3, label="LOWESS", zorder=2)
	ax.plot(step_list, efs_list, '-', c='blue', label="Raw Efs", zorder=1)

	plt.legend(loc=8, prop={'size': 24})
	# plt.show()
	file_save = "efs_"+file_print.replace("txt", "png")
	if (save):
		plt.savefig("figures/"+file_save, bbox_inches='tight', dpi=600)

	# Eff plot 
	fig, ax = plt.subplots()
	fig.set_size_inches(9, 9)

	ax.set_ylabel(r'Energy per particle (kT)', fontsize=30)
	ax.set_xlabel(r'Time (step)', fontsize=30)
	ax.tick_params(axis='both', which='major', labelsize=26)
	ax.tick_params(axis='both', which='minor', labelsize=20)
	ax.xaxis.offsetText.set_fontsize(22)
	plt.gcf().subplots_adjust(left=0.15)


	lowess = sm.nonparametric.lowess(eff_list, step_list, frac=.09)
	lowess_x = list(zip(*lowess))[0]
	lowess_y = list(zip(*lowess))[1]
	f = interp1d(lowess_x, lowess_y, bounds_error=False)
	xnew = np.linspace(step_list[0], step_list[-1], 10000)
	ynew = f(xnew)
	ax.plot(xnew, ynew, '-', c='red', linewidth=3, label="LOWESS", zorder=2)
	ax.plot(step_list, eff_list, '-', c='blue', label="Raw Eff", zorder=1)

	plt.legend(loc=4, prop={'size': 24})
	# plt.show()
	file_save = "eff_"+file_print.replace("txt", "png")
	if (save):
		plt.savefig("figures/"+file_save, bbox_inches='tight', dpi=600)


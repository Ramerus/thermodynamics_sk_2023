#!/usr/bin/env python3
from os import listdir
from os.path import isfile, join
import sys
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy.interpolate import interp1d
from matplotlib.transforms import blended_transform_factory
from scipy import interpolate
from scipy.optimize import curve_fit


def gamma_bet(value, gamma_m, k1_k2):

	value_cur = gamma_m * value*k1_k2/(1+(k1_k2-1)*value)/(1-value)

	return value_cur


if __name__ == "__main__":

	gamma_all = []	 	# gammas for the final plot 
	rho_all = []    	# rhos for the final plot 
	path = "./"
	size = 30			# size of a lattice (one side)
	rho_vle = 0.0185	# equilibrium density
	gt_crit = 1.0e7 	# only this values for averagine 
	save = 0 	# Tag for saving an image for files with positions 
	files = [f for f in listdir(path) if isfile(join(path, f))]
	energy = sys.argv[1] 
	files_position = []
	files_gamma = []
	for file in files:
		if ("positions_" in file) and (energy in file) :
			files_position.append(file)
			file_gamma = file.replace("positions", "gamma")
			files_gamma.append(file_gamma)
		
	if (len(files_position) == 0):
		raise ValueError("No files for selected energy")

	files_gamma.sort()
	files_position.sort()
	files_save = [1, 15, 29]
	files_save = []
	for idx_files, file in enumerate(files_position):
		rho_cur = file.split("_")
		rho_cur = rho_cur[-1]
		rho_cur = rho_cur.replace(".txt", "")
		rho_cur = float(rho_cur)
		rho_all.append(rho_cur/rho_vle)

		with open(file, "r") as f:
			lines = f.readlines()

		arr = np.zeros((size, size))
		x_coords = []
		y_coords = [] 
		for line in lines:
			line_c = line.split()
			x = int(line_c[0])
			y = int(line_c[1])
			arr[x][y] = 1

		if (idx_files in files_save):
			fig, ax = plt.subplots()
			ax.tick_params(axis='both', which='major', labelsize=16)
			ax.tick_params(axis='both', which='minor', labelsize=12)
			plt.imshow(arr)
			if (not save):
				plt.show()
			file_save = file.replace("txt", "png")
			if (save):
				plt.savefig("figures/"+file_save,bbox_inches='tight', dpi=150)


		# GAMMA FROM STEP
		with open(files_gamma[idx_files], "r") as f:
			lines = f.readlines()

		step_list = []
		gamma_list = []
		for line in lines:
			line_c = line.split()
			step_cur = int(line_c[1])
			gamma_cur = float(line_c[0])
			gamma_list.append(gamma_cur)
			step_list.append(step_cur)
		
		idx_sum = 0
		for idx, step in enumerate(step_list):
			if (step > gt_crit):
				idx_sum = idx - 1
				break
		gamma_mean = sum(gamma_list[idx_sum:]) / len(gamma_list[idx_sum:])
		gamma_all.append(gamma_mean)
		
		if (idx_files in files_save):
			fig, ax = plt.subplots()
			ax.tick_params(axis='both', which='major', labelsize=16)
			ax.tick_params(axis='both', which='minor', labelsize=12)
			ax.set_ylabel(r'$\Gamma$', fontsize=22)
			ax.set_ylabel(r'Step (number)', fontsize=22)
			plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
			plt.plot(step_list, gamma_list, '-', c='blue')
			ylim = ax.get_ylim()
			xlim = ax.get_xlim()
			plt.plot([gt_crit, gt_crit], ylim, '--', c='red', linewidth=3)

			text_annotate = r"Mean $\Gamma(\rho_{\mathrm{bulk}})$ = "
			text_annotate += str(round(gamma_mean,3))
			tform = blended_transform_factory(ax.transData, ax.transAxes)
			plt.text(xlim[1]*0.5, ylim[1]*0.05, text_annotate, fontsize=14)
			if (not save):
				plt.show()
			file_save = files_gamma[idx_files].replace("txt", "png")
			if (save):
				plt.savefig("figures/"+file_save,bbox_inches='tight', dpi=150)
		
		
# Final curve for the Gamma from rho_bulk 
gamma_all = [x for _,x in sorted(zip(rho_all,gamma_all))]
gamma_all_2 = [0.0, 0.0066666612500000015, 0.0070833300000000005, 0.01583332125, 0.031249988749999992, 
			0.0300000025, 0.04208331624999998, 0.048749989999999944, 0.050416677499999965,
			 0.07624997124999994, 0.09208329374999998, 0.10791666499999997, 0.13666664874999998, 
			 0.15208330000000006, 0.18958329625000006, 0.20000000625000003, 0.2295833587500001,
			  0.28916667125000006, 0.3287499637499999, 0.33541663000000005, 0.44708330875, 
			  0.46916668375000004, 0.5883333124999999, 0.6799999000000001, 0.75999998375, 
			  0.9908338250000004, 1.1966665000000003, 1.4333331249999999, 1.6924999124999995, 2.727500075]
rho_all.sort()

fig, ax = plt.subplots()
ax.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(axis='both', which='minor', labelsize=12)
ax.set_ylabel(r'$\Gamma$', fontsize=22)
ax.set_xlabel(r'$\rho$/$\rho_{\mathrm{bulk}}$ ', fontsize=22)
ax.scatter(rho_all, gamma_all, c='red', edgecolor='black', s=40, marker='o', zorder=1, label=r"GCMC $E_{\mathrm{sf}}=-2.0$ kT ")
# fitting with BET model
p_initial = [1, 5.0]
fit_l = 0
fit_r = 28
popt, _pcov = curve_fit(gamma_bet, rho_all[fit_l:fit_r], gamma_all[fit_l:fit_r], p0 = p_initial)
fit_ranges = np.linspace(rho_all[fit_l], rho_all[fit_r], 400)
gamma_m_4 = popt[0]
k1_k2_4 = popt[1]
print(popt)
# gamma_m_4 = 0.6
# k1_k2_4 = 15.0
curve_int = gamma_bet(fit_ranges, gamma_m_4, k1_k2_4)
ax.plot(fit_ranges, curve_int, '-', c='b', linewidth=2, label=r'BET $E_{\mathrm{sf}}=-4.0$ kT', zorder=3)

ax.scatter(rho_all, gamma_all_2, c='red', edgecolor='black', s=40, marker='s', zorder=1, label=r"GCMC $E_{\mathrm{sf}}=-4.0$ kT ")
# fitting with BET model
p_initial = [1, 0.1]
popt, _pcov = curve_fit(gamma_bet, rho_all, gamma_all_2, p0 = p_initial)
fit_ranges = np.linspace(rho_all[0], rho_all[-1], 400)
print(popt)
gamma_m_2 = popt[0]
k1_k2_2 = popt[1]
curve_int = gamma_bet(fit_ranges, gamma_m_2, k1_k2_2)
ax.plot(fit_ranges, curve_int, '-', c='g', linewidth=2, label=r'BET $E_{\mathrm{sf}}=-2.0$ kT', zorder=3)

ylim = ax.get_ylim()
xlim = ax.get_xlim()
plt.plot([1,1], ylim, '--', c='black')
plt.legend(loc=2, prop={'size': 12})

if (not save):
	plt.show()
file_save = "final_curve_fit_" + energy + ".png"
if (save):
	plt.savefig("figures/"+file_save,bbox_inches='tight', dpi=150)

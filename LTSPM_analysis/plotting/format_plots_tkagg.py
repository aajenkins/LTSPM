# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 10:55:17 2017

@author: alec
"""


def format_plots(plt, square=False, padding=0, small=True, tight=True):
	figNums = plt.get_fignums()
	flen = len(figNums)

	if (small):
		fheight = int((900-(3*24))/2)
		if (square):
			fwidth = fheight
		else:
			fwidth = int(1440/3) - padding
		for i in range(flen):
			plt.figure(figNums[i])
			fwin = plt.get_current_fig_manager().window
			x =  (padding + fwidth) * (i%3)
			y = 22 + (23 + padding + fheight) * (int(i/3)%2)
			fwin.geometry('%dx%d+%d+%d'%(fwidth, fheight, x, y))
			if (tight):
				plt.gcf().set_tight_layout(True)
			# print(fwidth, fheight, x, y)

	else:
		fheight = 650
		fwidth = 900
		for i in range(flen):
			plt.figure(figNums[i])
			fwin = plt.get_current_fig_manager().window
			x =  100 * (i%8)
			y = 22
			fwin.geometry('%dx%d+%d+%d'%(fwidth, fheight, x, y))
			if (tight):
				plt.gcf().set_tight_layout(True)

def format_plot(plt, fwidth=700, fheight=600, tight=True, x =  10, y = 22):
	fwin = plt.get_current_fig_manager().window
	fwin.geometry('%dx%d+%d+%d'%(fwidth, fheight, x, y))
	if (tight):
		plt.gcf().set_tight_layout(True)

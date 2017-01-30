# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 10:55:17 2017

@author: alec
"""

def format_plot(plt, width, height, wx, wy, tight = True):
	# fig = plt.figure()
	if tight:
		plt.gcf().set_tight_layout(True)
	plt.get_current_fig_manager().window.setGeometry(wx, wy, width, height)
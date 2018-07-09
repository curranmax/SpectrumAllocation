
import argparse
import matplotlib.pyplot as plt
from numpy import median
import math

from runner import readInParamResult, runExperiment

class Point:
	def __init__(self, x, y, group):
		self.x = x
		self.y = y
		self.group = group

def groupValues(raw_values, xv, xfunc, yv, yfunc, gv, fv):
	points = []
	for param, result in raw_values:
		skip = False
		for f_var in fv:
			this_v = param.getValue(f_var)
			if this_v not in fv[f_var]:
				skip = True

		if skip:
			continue

		this_x = xfunc(*tuple(getattr(param, x) for x in xv))

		this_group = tuple(getattr(param, g) for g in gv)

		tmp_y = tuple((result.getValue(y) if result.hasAttr(y) else param.getValue(y)) for y in yv)
		if any(isinstance(k, list) for k in tmp_y):
			lens = [len(k) for k in tmp_y if isinstance(k, list)]
			if any(v != lens[0] for v in lens):
				raise Exception('All lists must be of same length')

			tmp_y = [(k if isinstance(k, list) else [k] * lens[0]) for k in tmp_y]

			for vs in zip(*tmp_y):
				if None not in vs:
					this_y = yfunc(*vs)
					points.append(Point(this_x, this_y, this_group))
		else:
			if None not in tmp_y:
				this_y = yfunc(*tmp_y)
				points.append(Point(this_x, this_y, this_group))

	groups = set(p.group for p in points)
	x_by_group = {group:set() for group in groups}
	for p in points:
		x_by_group[p.group].add(p.x)

	data = {group:{x:[] for x in x_by_group[group]} for group in groups}
	for p in points:
		data[p.group][p.x].append(p.y)
	return data

def plotData(data, reduce_function = None, include_errbars = False, scatter_plot = False, title = None, xlabel = None, ylabel = None):
	handles = []
	legend_labels = []
	for group in data:
		_xvs = [x for x in sorted(data[group])]
		_yvs = [reduce_function(data[group][xv]) for xv in _xvs]

		xvs = []
		yvs = []
		for xv, yv in zip(_xvs, _yvs):
			if isinstance(yv, list):
				for yv_element in sorted(yv):
					xvs.append(xv)
					yvs.append(yv_element)
			else:
				xvs.append(xv)
				yvs.append(yv)

		if scatter_plot:
			h = plt.scatter(xvs, yvs, alpha = 0.25)
		elif include_errbars:
			yerr_low = [yv - min(data[group][xv]) for xv, yv in zip(xvs, yvs)]
			yerr_high = [max(data[group][xv]) - yv for xv, yv in zip(xvs, yvs)]

			h, _, _ = plt.errorbar(xvs, yvs, yerr = [yerr_low, yerr_high])
		else:
			h, = plt.plot(xvs, yvs)

		handles.append(h)
		legend_labels.append(str(group))

	if title != None:
		plt.title(title)

	if xlabel != None:
		plt.xlabel(xlabel)

	if ylabel != None:
		plt.ylabel(ylabel)

	if len(data) > 1:
		plt.legend(handles, legend_labels)

	plt.show()

pos_colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', '0.5']
def getColor(k, colors):
	if k in colors:
		return colors[k]
	for c in pos_colors:
		if all(colors[x] != c for x in colors):
			colors[k] = c
			return c
	raise Exception('NO MORE COLORS LEFT!!')

def makeBoxPlot(data, reduce_function = None, include_errbars = False, scatter_plot = False, title = None, xlabel = None, ylabel = None):
	handles = []
	legend_label = []

	fig = plt.figure()
	ax = plt.axes()

	# Figures Out x-coords for the inidividual box-plots
	# intra_buffer and inter_buffer are constants that determine the layout
	intra_buffer = .1
	inter_buffer = .3
	num_plots = len(data)
	plot_width = (1.0 - inter_buffer - intra_buffer * (num_plots - 1)) / num_plots
	positions = {k:[] for k in data}
	include_fliers = True
	for i, k in enumerate(sorted(data)):
		start_position = .5 + inter_buffer / 2.0 + plot_width / 2.0 + i * (plot_width + intra_buffer)
		positions[k] = [start_position + i for i in range(len(data[k]))]

	x_labels = {}
	colors = {}
	for k in sorted(data):
		x_vals = sorted([xv for xv in data[k] if len(data[k][xv]) > 0])
		for i, xv in enumerate(x_vals):
			x_labels[(i + 1)] = xv
		y_vals = [[yv for yv in data[k][xv]] for xv in x_vals]
		try:
			bp = plt.boxplot(y_vals, positions = positions[k], widths = plot_width, sym = ('+' if include_fliers else ''))
		except ValueError:
			bp = None

		this_color = getColor(k, colors)
		if bp != None:
			for elem in ('boxes','caps','whiskers', 'fliers', 'medians'):
				[plt.setp(bp[elem][idx], color = this_color) for idx in range(len(bp[elem]))]

	plt.xlabel(xlabel)
	plt.xlim([.5 - inter_buffer / 2.0, len(x_labels) + .5 + inter_buffer / 2.0])
	ax.set_xticklabels([x_labels[i] for i in sorted(x_labels)])
	ax.set_xticks(sorted(x_labels.keys()))
	plt.ylabel(ylabel)
	# plt.ylim([-1.0, 5])
	plt.title(title)

	# Make legend]
	if len(data) > 1:
		for k in sorted(data):
			h, = plt.plot([1,1], color = colors[k], linestyle = '-')
			handles.append(h)
			legend_label.append(k)
		plt.legend(handles, legend_label)
		for h in handles:
			h.set_visible(False)

	plt.show()

def tryConvert(val):
	if val == 'True':
		return True
	if val == 'False':
		return False

	try:
		return int(val)
	except ValueError:
		pass

	try:
		return float(val)
	except ValueError:
		pass

	return val

def percentDifference(v1, v2, unit_type):
	if unit_type == 'db':
		v1 = pow(10.0, v1 / 10.0)
		v2 = pow(10.0, v2 / 10.0)
	return float(abs(v1 - v2)) / float(abs(v2)) * 100.0

def dbDifference(v1, v2, unit_type):
	if unit_type == 'abs':
		if v1 <= 0.0:
			v1 = math.pow(10.0, -50.0)
		v1 = 10.0 * math.log(v1, 10.0)
		v2 = 10.0 * math.log(v2, 10.0)
	return abs(v1 - v2)

def average(vals):
	return sum(vals) / float(len(vals))

def allPoints(vals):
	return vals

x_labels = {('num_ss_selection',): 'Number of SS Selected',
			('ss_receive_power_alpha',): 'Receive Power Alpha',
			('ss_path_loss_alpha',): 'Path Loss Alpha',
			('num_ss',): 'Total Number of SS',
			('num_pu',): 'Total Number of PU'}

def getXLabel(xv):
	if xv in x_labels:
		return x_labels[xv]
	else:
		print 'Unrecognized xv:', xv
		return str(xv)

y_labels = {'percent_diff_secure_vs_plain': 'Percent Difference between Secure and Plain',
			'percent_diff_plain_vs_ground': 'Percent Difference between Plain and Ground',
			'percent_diff_secure_vs_ground': 'Percent Difference between Secure and Ground',
			'db_diff_secure_vs_plain': 'dB Difference between Secure and Plain',
			'db_diff_plain_vs_ground': 'dB Difference between Plain and Ground',
			'db_diff_secure_vs_ground': 'dB Difference between Secure and Ground',
			'time_per_request': 'Time per SU Request (seconds)',
			'preprocess_time': 'Pre-process time (seconds)'}

def getYLabel(yv):
	if yv in y_labels:
		return y_labels[yv]
	else:
		print 'Unrecognized yv:', yv
		return str(yv)

rfuncs = {'average': average, 'min': min, 'max': max, 'median': median, 'all': allPoints}

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Graphs output of S2-PC experiments.')

	# X values
	x_values = [('num_pu', 'npu'), ('num_su', 'nsu'), ('num_ss', 'nss'),
				('location_range', 'lr'), ('ld_path_loss0', 'ld_pl0'), ('ld_dist0', 'ld_d0'),
				('num_ss_selection', 'nss_s'), ('num_pu_selection', 'npu_s'), ('ss_receive_power_alpha', 'rpa'), ('ss_path_loss_alpha', 'pla'),
				('num_float_bits', 'nfb'), ('s2_pc_bit_count', 'bc'),
				('grid_x', 'gx'), ('grid_y', 'gy')]
	for full_xv, short_xv in x_values:
		parser.add_argument('-x_' + short_xv, '--use_' + full_xv + '_for_x_value', action = 'store_true', help = 'If given, uses ' + full_xv + ' for the x dimension.')

	# Y values
	y_values = [('preprocess_time' ,'ppt'), ('time_per_request', 'tpr'),
				('percent_diff_secure_vs_plain', 'svp'), ('percent_diff_plain_vs_ground', 'pvg'), ('percent_diff_secure_vs_ground', 'svg'),
				('db_diff_secure_vs_plain', 'svp_db'), ('db_diff_plain_vs_ground', 'pvg_db'), ('db_diff_secure_vs_ground', 'svg_db')]
	for full_yv, short_yv in y_values:
		parser.add_argument('-y_' + short_yv, '--use_' + full_yv + '_for_y_value', action = 'store_true', help = 'If given, uses ' + full_yv + ' for the y dimension.')

	# Group values
	g_values = x_values + [('unit_type', 'ut'), ('algo_order', 'ao'), ('selection_algo', 'sa')]
	for full_gv, short_gv in g_values:
		parser.add_argument('-g_' + short_gv, '--use_' + full_gv + '_for_group_value', action = 'store_true', help = 'If given, uses ' + full_gv + ' for grouping data.')

	# Filter values
	f_values = list(g_values)
	for full_fv, short_fv in f_values:
		parser.add_argument('-f_' + short_fv, '--use_' + full_fv + '_for_filter_value', metavar = 'FILTER_VALUE', type = str, nargs = '+', default = None , help = 'Only plots points with values of ' + full_gv + ' that match one of the given values.')

	parser.add_argument('-reduce', '--reduce_function', metavar = 'FUNCTION', type = str, nargs = 1, default = ['average'], help = 'Function to reduce multiple functions. Can be (' + ', '.join(rfuncs.keys()) + ')')
	parser.add_argument('-errbar', '--include_errbars', action = 'store_true', help = 'If given, includes error bars in the graph.')
	parser.add_argument('-scatter', '--make_scatter_plot', action = 'store_true', help = 'If given, creates a scatter plot instead of line plot.')
	parser.add_argument('-box_plot', '--make_box_plot', action = 'store_true', help = 'If given, creates a box plot instead of line plot.')

	parser.add_argument('-t', '--title', metavar = 'TITLE', type = str, nargs = 1, default = [''], help = 'Title of the plot')

	parser.add_argument('-in', '--in_files', metavar = 'IN_FILE', type = str, nargs = '+', default = [], help = 'List of files to get data from.')

	args = parser.parse_args()

	raw_values = readInParamResult(args.in_files)

	xv = None
	xfunc = None
	for full_xv, _ in x_values:
		if getattr(args, 'use_' + full_xv + '_for_x_value'):
			if xv != None:
				raise Exception('Cannot specify multiply x values')
			xv = (full_xv,)
			xfunc = lambda x: x

	yv = None
	yfunc = None
	complex_y_values = {'percent_diff_secure_vs_plain': (('su_transmit_power', 'secure'), ('su_transmit_power', 'plain'), 'unit_type'),
						'percent_diff_plain_vs_ground': (('su_transmit_power', 'plain'), ('su_transmit_power', 'ground'), 'unit_type'),
						'percent_diff_secure_vs_ground': (('su_transmit_power', 'secure'), ('su_transmit_power', 'ground'), 'unit_type'),
						'db_diff_secure_vs_plain': (('su_transmit_power', 'secure'), ('su_transmit_power', 'plain'), 'unit_type'),
						'db_diff_plain_vs_ground': (('su_transmit_power', 'plain'), ('su_transmit_power', 'ground'), 'unit_type'),
						'db_diff_secure_vs_ground': (('su_transmit_power', 'secure'), ('su_transmit_power', 'ground'), 'unit_type')}
	complex_y_funcs = {'percent_diff_secure_vs_plain': percentDifference,
						'percent_diff_plain_vs_ground': percentDifference,
						'percent_diff_secure_vs_ground': percentDifference,
						'db_diff_secure_vs_plain': dbDifference,
						'db_diff_plain_vs_ground': dbDifference,
						'db_diff_secure_vs_ground': dbDifference}
	this_y_value = None
	for full_yv, _ in y_values:
		if getattr(args, 'use_' + full_yv + '_for_y_value'):
			if yv != None:
				raise Exception('Cannot specify multiply y values')
			if full_yv in complex_y_values:
				yv = complex_y_values[full_yv]
				yfunc = complex_y_funcs[full_yv]
			else:
				yv = (full_yv,)
				yfunc = lambda y: y
			this_y_value = full_yv

	gv = []
	for full_gv, _ in g_values:
		if getattr(args, 'use_' + full_gv + '_for_group_value'):
			gv.append(full_gv)

	fv = {}
	for full_fv, _ in f_values:
		vals = getattr(args, 'use_' + full_fv + '_for_filter_value')
		if vals != None and len(vals) > 0:
			fv[full_fv] = map(tryConvert, vals)

	data = groupValues(raw_values, xv, xfunc, yv, yfunc, gv, fv)

	if args.make_box_plot:
		makeBoxPlot(data, title = args.title[0], xlabel = getXLabel(xv), ylabel = getYLabel(this_y_value))
	else:
		plotData(data, reduce_function = rfuncs[args.reduce_function[0]], include_errbars = args.include_errbars, scatter_plot = args.make_scatter_plot,
				title = args.title[0], xlabel = getXLabel(xv), ylabel = getYLabel(this_y_value))

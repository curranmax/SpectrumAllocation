
import argparse
import matplotlib.pyplot as plt
from numpy import median

from runner import readInParamResult, runExperiment

class Point:
	def __init__(self, x, y, group):
		self.x = x
		self.y = y
		self.group = group

def groupValues(raw_values, xv, xfunc, yv, yfunc, gv):
	points = []
	for param, result in raw_values:
		this_x = xfunc(*tuple(getattr(param, x) for x in xv))

		this_group = tuple(getattr(param, g) for g in gv)

		tmp_y = tuple((result.getValue(y) if result.hasAttr(y) else param.getValue(y)) for y in yv)
		if any(isinstance(k, list) for k in tmp_y):
			lens = [len(k) for k in tmp_y if isinstance(k, list)]
			if any(v != lens[0] for v in lens):
				raise Exception('All lists must be of same length')

			tmp_y = [(k if isinstance(k, list) else [k] * lens[0]) for k in tmp_y]

			for vs in zip(*tmp_y):
				this_y = yfunc(*vs)
				points.append(Point(this_x, this_y, this_group))
		else:
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

def plotData(data, reduce_function = None, include_errbars = False, title = None, xlabel = None, ylabel = None):

	handles = []
	legend_labels = []
	for group in data:
		xvs = [x for x in sorted(data[group])]
		yvs = [reduce_function(data[group][xv]) for xv in xvs]

		if include_errbars:
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

def percentDifference(v1, v2, unit_type):
	if unit_type == 'db':
		v1 = pow(10.0, v1 / 10.0)
		v2 = pow(10.0, v2 / 10.0)
	return float(abs(v1 - v2)) / float(abs(v2)) * 100.0

def dbDifference(v1, v2, unit_type):
	if unit_type == 'abs':
		v1 = 10.0 * log(v1, 10.0)
		v2 = 10.0 * log(v2, 10.0)
	return abs(v1 - v2)

def average(vals):
	return sum(vals) / float(len(vals))

x_labels = {('num_ss_selection',): 'Number of SS Selected',
			('ss_receive_power_alpha',): 'Receive Power Alpha',
			('ss_path_loss_alpha',): 'Path Loss Alpha'}

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

rfuncs = {'average': average, 'min': min, 'max': max, 'median': median}

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Graphs output of S2-PC experiments.')

	# X values
	x_values = [('num_pu', 'npu'), ('num_su', 'nsu'), ('num_ss', 'nss'),
				('location_range', 'lr'), ('ld_path_loss0', 'ld_pl0'), ('ld_dist0', 'ld_d0'),
				('num_ss_selection', 'nss_s'), ('ss_receive_power_alpha', 'rpa'), ('ss_path_loss_alpha', 'pla'),
				('num_float_bits', 'nfb'), ('s2_pc_bit_count', 'bc')]
	for full_xv, short_xv in x_values:
		parser.add_argument('-x_' + short_xv, '--use_' + full_xv + '_for_x_value', action = 'store_true', help = 'If given, uses ' + full_xv + ' for the x dimension.')

	# Y values
	y_values = [('preprocess_time' ,'ppt'), ('time_per_request', 'tpr'),
				('percent_diff_secure_vs_plain', 'svp'), ('percent_diff_plain_vs_ground', 'pvg'), ('percent_diff_secure_vs_ground', 'svg'),
				('db_diff_secure_vs_plain', 'svp_db'), ('db_diff_plain_vs_ground', 'pvg_db'), ('db_diff_secure_vs_ground', 'svg_db')]
	for full_yv, short_yv in y_values:
		parser.add_argument('-y_' + short_yv, '--use_' + full_yv + '_for_y_value', action = 'store_true', help = 'If given, uses ' + full_yv + ' for the y dimension.')

	# Group values
	g_values = x_values + []
	for full_gv, short_gv in g_values:
		parser.add_argument('-g_' + short_gv, '--use_' + full_gv + '_for_group_value', action = 'store_true', help = 'If given, uses ' + full_gv + ' for grouping data.')

	parser.add_argument('-reduce', '--reduce_function', metavar = 'FUNCTION', type = str, nargs = 1, default = ['average'], help = 'Function to reduce multiple functions. Can be (' + ', '.join(rfuncs.keys()) + ')')
	parser.add_argument('-errbar', '--include_errbars', action = 'store_true', help = 'If given, includes error bars in the graph.')

	parser.add_argument('-t', '--title', metavar = 'TITLE', type = str, nargs = 1, default = [''], help = 'Title of the plot')

	parser.add_argument('-in', '--in_files', metavar = 'IN_FILE', type = str, nargs = '+', default = [], help = 'List of files to get data from.')

	args = parser.parse_args()

	raw_values = readInParamResult(args.in_files)

# 	vals = [(percentDifference(plain, ground), param, result) for param, result in raw_values
# 				for plain, ground in zip(result.su_transmit_power['plain'], result.su_transmit_power['ground'])]
# 	v, param, result = max(vals)

# 	param.rand_seed = result.rand_seed

# 	runExperiment(param, no_run = True)
# else:

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

	data = groupValues(raw_values, xv, xfunc, yv, yfunc, gv)
	plotData(data, reduce_function = rfuncs[args.reduce_function[0]], include_errbars = args.include_errbars,
				title = args.title[0], xlabel = getXLabel(xv), ylabel = getYLabel(this_y_value))

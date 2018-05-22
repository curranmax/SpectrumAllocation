
import subprocess
import argparse
from datetime import datetime
import itertools
from operator import mul
from copy import deepcopy

# Experiments
TEST = 'test'

SMALL_LD_VARY_NUM_SS_SELECT = 'vary_ss_select'
SMALL_LD_VARY_RP_ALPHA = 'vary_rp_alpha'
SMALL_LD_VARY_PL_ALPHA = 'vary_pl_alpha'
SMALL_LD_VARY_RP_AND_PL_ALPHA = 'vary_rp_and_pl_alpha'

SMALL_LD_VARY_NUM_SS_SELECT_AND_ALPHAS = 'vary_ss_select_and_alphas'

ALL_EXPERIMENT_IDS = [TEST,
						SMALL_LD_VARY_NUM_SS_SELECT,
						SMALL_LD_VARY_RP_ALPHA, SMALL_LD_VARY_PL_ALPHA, SMALL_LD_VARY_RP_AND_PL_ALPHA,
						SMALL_LD_VARY_NUM_SS_SELECT_AND_ALPHAS]

# Change parameters
NUM_SS_SELECTION = 'num_ss_selection'
RP_ALPHA = 'ss_receive_power_alpha'
PL_ALPHA = 'ss_path_loss_alpha'

def makeDefaultExperimentParam():
	param = ExperimentParam()

	param.num_pu = 5
	param.num_su = 10
	param.num_ss = 100

	param.location_range = 100

	param.propagation_model = 'log_distance'

	param.ld_path_loss0 = 1
	param.ld_dist0 = 5
	param.ld_gamma = 1

	param.num_ss_selection = 0
	algo_order = 'split_then_idw'
	path_loss_type = 'ratio'
	param.ss_receive_power_alpha = 1
	param.ss_path_loss_alpha = 1

	param.num_float_bits = 16
	param.s2_pc_bit_count = 64

	return param

class ExperimentParam:
	def __init__(self):
		self.rand_seed = None

		self.num_pu = None
		self.num_su = None
		self.num_ss = None

		self.location_range = None

		self.propagation_model = None

		self.ld_path_loss0 = None
		self.ld_dist0 = None
		self.ld_gamma = None

		self.num_ss_selection = None

		self.algo_order = None
		self.path_loss_type = None

		self.ss_receive_power_alpha = None
		self.ss_path_loss_alpha = None

		self.num_float_bits = None
		self.s2_pc_bit_count = None

def makeExperimentResult(vals):
	result = ExperimentResult()
	for val in vals:
		if '|' not in val:
			print 'Could\'t parse val:', val
			continue

		var_name, val_type, val = val.split('|')

		var_name, val = convert(var_name, val_type, val)
		vars(result)[var_name] = val

	return result

class ExperimentResult:
	def __init__(self):
		preprocess_time = None
		su_transmit_power = None
		rand_seed = None
		time_per_request = None

	def getCols(self):
		# each value is a 3 tuple (var_name, var_name_for_output, val_type)
		cols = []
		for k in vars(self):
			if vars(self)[k] == None:
				continue

			if isinstance(getattr(self, k), int):
				x = (k, None, 'int')
			elif isinstance(getattr(self, k), float):
				x = (k, None, 'float')
			elif isinstance(getattr(self, k), str):
				x = (k, None, 'str')
			elif isinstance(getattr(self, k), dict) \
					and all(isinstance(label, str) and isinstance(vals, list) \
							and all(isinstance(v, float) for v in vals) \
						for label, vals in getattr(self, k).iteritems()):
				x = (k, tuple(label for label in getattr(self, k)), 'list(' + ','.join(['float' for label in getattr(self, k)]) + ')')
			else:
				raise Exception('Can\'t determine the type of: ' + str(getattr(self, k)))
			cols.append(x)
		return cols

	def getStr(self, var_name, var_order):
		if var_order == None:
			return str(getattr(self, var_name))
		
		return ','.join(':'.join(map(str, vs)) for vs in zip(*tuple(getattr(self, var_name)[label] for label in var_order)))

	def getValue(self, yv):
		if isinstance(yv, str):
			return getattr(self, yv)
		if isinstance(yv, tuple) and len(yv) == 2:
			return getattr(self, yv[0])[yv[1]]

		raise Exception('Invalid attribute: ' + str(yv))

def writeOutParamResult(filename, out_values):
	f = open(filename, 'a')
	
	shared_param, shared_param_str = getSharedParam(out_values)
	f.write('PARAM ' + shared_param_str + '\n')

	param_cols = [(k, getType(getattr(shared_param, k)[0])) for k in vars(shared_param) if isinstance(getattr(shared_param, k), list)]
	result_cols = list(set(col for _, result in out_values for col in result.getCols()))
	f.write('COL_NAMES PARAM_COL_NAMES ' + ' '.join('|'.join(v) for v in param_cols) + \
			' RESULT_COL_NAMES ' + ' '.join('|'.join((var_name + ('' if var_order == None else '(' + ','.join(var_order) + ')'), val_type)) for var_name, var_order, val_type in result_cols) + '\n')

	for param, result in out_values:
		line = []
		for k, _ in param_cols:
			line.append(str(getattr(param, k)))

		for var_name, var_order, _ in result_cols:
			line.append(result.getStr(var_name, var_order))

		f.write(' '.join(line) + '\n')

def readInParamResult(filenames):
	rv = []

	for filename in filenames:
		f = open(filename, 'r')

		cur_params = None
		cur_param_cols = None
		cur_result_cols = None

		for line in f:
			spl = line.split()
			if len(spl) == 0:
				continue

			if spl[0] == 'PARAM':
				cur_params = ExperimentParam()
				for v in spl[1:]:
					var_name, val_type, val_str = v.split('|')
					var_name, val = convert(var_name, val_type, val_str)
					vars(cur_params)[var_name] = val
			elif spl[0] == 'COL_NAMES':

				cur_param_cols  = [tuple(v.split('|')) for v in spl[spl.index('PARAM_COL_NAMES') + 1 : spl.index('RESULT_COL_NAMES')]]
				cur_result_cols = [tuple(v.split('|')) for v in spl[spl.index('RESULT_COL_NAMES') + 1:]]
			else:
				if cur_params == None or cur_param_cols == None or cur_result_cols == None:
					raise Exception('Invalid file syntax')

				this_params = deepcopy(cur_params)
				this_result = ExperimentResult()

				for (var_name, val_type), val_str in zip(cur_param_cols, spl[:len(cur_param_cols)]):
					var_name, val = convert(var_name, val_type, val_str)
					vars(this_params)[var_name] = val

				for (var_name, val_type), val_str in zip(cur_result_cols, spl[len(cur_param_cols):]):
					var_name, val = convert(var_name, val_type, val_str)
					vars(this_result)[var_name] = val

				rv.append((this_params, this_result))
	return rv

def runExperiment(param):
	args = ['../s2pc', '-brief_out']
	vs =[(param.rand_seed, 'rand_seed'),
			(param.num_pu, 'npu'), (param.num_ss, 'nss'), (param.num_su, 'nsu'), (param.location_range, 'lr'),
			(param.propagation_model, 'pm'),
			(param.ld_path_loss0, 'ld_pl0'), (param.ld_dist0, 'ld_d0'), (param.ld_gamma, 'ld_g'),
			(param.num_ss_selection, 'nss_s'),
			(param.algo_order, 'ao'), (param.path_loss_type, 'plt'),
			(param.ss_receive_power_alpha, 'rpa'), (param.ss_path_loss_alpha, 'pla'),
			(param.num_float_bits, 'float_bits'), (param.s2_pc_bit_count, 'bit_count')]

	for val, flag in vs:
		if val != None:
			args += ['-' + flag, str(val)]

	rv = None
	for i in range(5):
		process = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		process.wait()

		# If there is any output to stderr, its assumed that something went wrong, and the experiment must be run again.
		err_out = [v for v in process.stderr]
		if len(err_out) != 0:
			continue

		rv = makeExperimentResult([v[:-1] for v in process.stdout])
		break

	if rv == None:
		raise Exception('Unable to run program, got:\n' + ''.join([v for v in process.stderr]))

	return rv

def getSharedParam(out_values):
	shared_param = ExperimentParam()
	for param, _ in out_values:
		for k in vars(param):
			if callable(getattr(param, k)) or k.startswith("__") or getattr(param, k) == None:
				continue
			if getattr(shared_param, k) == None:
				vars(shared_param)[k] = getattr(param, k)
			elif isinstance(getattr(shared_param, k), list):
				if getattr(param, k) not in getattr(shared_param, k):
					vars(shared_param)[k].append(getattr(param, k))
			elif getattr(shared_param, k) != getattr(param, k):
				vars(shared_param)[k] = [getattr(shared_param, k), getattr(param, k)]

	return shared_param, ' '.join([paramToString(k, getattr(shared_param, k)) for k in vars(shared_param) if not(callable(getattr(param, k)) or k.startswith("__"))])

def paramToString(var_name, value):
	type_str = ''
	value_str = ''
	if value == None:
		type_str = 'None'
		value_str = 'None'
	elif isinstance(value, int):
		type_str = 'int'
		value_str = str(value)
	elif isinstance(value, float):
		type_str = 'float'
		value_str = str(value)
	elif isinstance(value, str):
		type_str = 'str'
		value_str = str(value)
	elif isinstance(value, list):
		if isinstance(value[0], int):
			type_str = 'list(int)'
			value_str = ','.join(map(str, value))
		elif isinstance(value[0], float):
			type_str = 'list(float)'
			value_str = ','.join(map(str, value))
		elif isinstance(value[0], str):
			type_str = 'list(str)'
			value_str = ','.join(map(str, value))
		else:
			raise Exception('Can\'t determine the type of: (' + str(var_name) + ', ' + str(value) + ')')

	return var_name + '|' + type_str + '|' + value_str

def getType(value):
	if value == None:
		return 'None'
	if isinstance(value, int):
		return 'int'
	if isinstance(value, float):
		return 'float'
	if isinstance(value, str):
		return 'str'
	raise Exception('Can\'t determine the type of: ' + str(value))

def convert(var_name, val_type, val_str):
	if val_type == 'None':
		return var_name, None
	if val_type == 'int':
		return var_name, int(val_str)
	if val_type == 'float':
		return var_name, float(val_str)
	if val_type == 'str':
		return var_name, val_str
	if val_type == 'list(int)':
		return var_name, list(map(int, val_str.split(',')))
	if val_type == 'list(float)':
		return var_name, list(map(int, val_str.split(',')))
	if val_type == 'list(str)':
		return var_name, list(map(str, val_str.split(',')))
	if val_type == 'list(float,float,float)':
		simple_var_name = var_name[:var_name.find('(')]
		val_labels = var_name[var_name.find('(') + 1 : var_name.find(')')].split(',')

		vals = map(lambda x: tuple(map(float, x.split(':'))), val_str.split(','))
		vals = {label:[v[i] for v in vals] for i, label in enumerate(val_labels)}
		return simple_var_name, vals

	raise Exception('Unknown type: (' + str(var_name) + ', ' + str(val_type) + ', ' + str(val_str) + ')')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Runs experiments. Has the parameters for each experiment saved, so they don\'t have to be reentered')
	
	current_time = datetime.now()

	parser.add_argument('-exp', '--experiment_identifier', metavar = 'EXPERIMENT_ID', type = str, nargs = '+', default = [], help = 'List of experiments to run. Must be value in (' + ', '.join(ALL_EXPERIMENT_IDS) + ')')
	parser.add_argument('-out', '--out_file', metavar = 'OUT_FILE', type = str, nargs = 1, default = [None], help = 'File to write data to')
	parser.add_argument('-nt', '--num_tests', metavar = 'NUM_TESTS', type = int, nargs = 1, default = [1], help = 'Number of tests to run')

	args = parser.parse_args()
	experiments = args.experiment_identifier

	changes = []

	for experiment in experiments:
		if experiment == TEST:
			changes.append({NUM_SS_SELECTION: [1], 'path_loss_type':['db']})
		if experiment == SMALL_LD_VARY_NUM_SS_SELECT:
			changes.append({NUM_SS_SELECTION: [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], 'num_pu': [1], PL_ALPHA: [2], RP_ALPHA: [2],
							'location_range': [250.0], 'num_ss': [625]})
		if experiment == SMALL_LD_VARY_RP_ALPHA:
			changes.append({RP_ALPHA: [1, 2, 3, 4]})
		if experiment == SMALL_LD_VARY_PL_ALPHA:
			changes.append({PL_ALPHA: [1, 2, 3, 4]})
		if experiment == SMALL_LD_VARY_RP_AND_PL_ALPHA:
			changes.append({RP_ALPHA: [1, 2, 3, 4], PL_ALPHA: [1, 2, 3, 4]})
		if experiment == SMALL_LD_VARY_NUM_SS_SELECT_AND_ALPHAS:
			changes.append({NUM_SS_SELECTION: [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
							RP_ALPHA: [1, 2, 3, 4], PL_ALPHA: [1, 2, 3, 4]})

	num_experiments = sum(reduce(mul, [len(vals) for _, vals in change.iteritems()]) for change in changes) * args.num_tests[0]

	out_values = []
	exp_num = 1
	for i in range(args.num_tests[0]):
		for change in changes:
			ns = change.keys()
			for vs in itertools.product(*tuple(vs for k, vs in change.iteritems())):
				print 'Running experiment', exp_num, 'of', num_experiments, 'with parameters:', ns, vs
				exp_num += 1

				param = makeDefaultExperimentParam()

				for n, v in zip(ns, vs):
					vars(param)[n] = v

				result = runExperiment(param)
				out_values.append((param, result))


	out_file = args.out_file[0]
	if out_file == None:
		out_file = current_time.strftime('data/out_%b-%d_%H:%M:%S_' + '-'.join(experiments) + '.txt')

	writeOutParamResult(out_file, out_values)
	
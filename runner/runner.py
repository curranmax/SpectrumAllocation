
import subprocess
import argparse
from datetime import datetime, timedelta
import itertools
from operator import mul
from copy import deepcopy
import math

# Experiments
TEST = 'test'

VARY_NUM_SS_SELECT = 'vary_ss_select'

TIMING_TEST = 'timing_test'
PATH_LOSS_TEST = 'path_loss_test'

MED_TEST = 'med_test'
LARGE_TEST = 'large_test'

TEST_OUTPUT= 'test_output'
OUTPUT_RUN = 'output'

ALL_EXPERIMENT_IDS = [TEST,
						VARY_NUM_SS_SELECT,
						TIMING_TEST, PATH_LOSS_TEST,
						MED_TEST, LARGE_TEST,
						TEST_OUTPUT, OUTPUT_RUN]

# Change parameters
NUM_SS_SELECTION = 'num_ss_selection'
RP_ALPHA = 'ss_receive_power_alpha'
PL_ALPHA = 'ss_path_loss_alpha'

def makeDefaultExperimentParam():
	param = ExperimentParam()

	param.skip_s2pc = False

	param.central_entities = 'two_sms'

	param.num_pu = 5
	param.num_su = 10
	param.num_ss = 100

	param.num_pr_per_pu = 1
	param.pr_range = 0.0

	param.location_range = 100

	param.unit_type = 'abs'

	param.propagation_model = 'longley_rice'

	param.ld_path_loss0 = None
	param.ld_dist0 = None
	param.ld_gamma = None

	param.splat_cmd = 'splat'

	param.ref_lat = 40.75
	param.ref_long = 73.25

	param.splat_dir = "../splat/"
	param.sdf_dir = "sdf/"
	param.return_dir = "../runner/"

	param.num_ss_selection = 0
	param.num_pu_selection = 0

	param.do_plaintext_split = True
	param.no_pr_thresh_update = False

	param.selection_algo = 'sort'
	param.secure_write_algo = 'proposed'

	param.grid_x = 0
	param.grid_y = 0

	param.ss_receive_power_alpha = 1
	param.ss_path_loss_alpha = 1

	param.num_float_bits = 16
	param.s2_pc_bit_count = 64

	return param

class ExperimentParam:
	def __init__(self):
		self.rand_seed = None
		self.skip_s2pc = None

		self.central_entities = None

		self.num_pu = None
		self.num_su = None
		self.num_ss = None

		self.num_pr_per_pu = None
		self.pr_range = None

		self.location_range = None

		self.out_filename = None

		self.unit_type = None

		self.propagation_model = None

		self.ld_path_loss0 = None
		self.ld_dist0 = None
		self.ld_gamma = None

		self.splat_cmd = None

		self.ref_lat = None
		self.ref_long = None

		self.splat_dir = None
		self.sdf_dir = None
		self.return_dir = None

		self.num_ss_selection = None
		self.num_pu_selection = None

		self.do_plaintext_split = None
		self.no_pr_thresh_update = None

		self.selection_algo = None
		self.secure_write_algo = None

		self.grid_x = None
		self.grid_y = None

		self.ss_receive_power_alpha = None
		self.ss_path_loss_alpha = None

		self.num_float_bits = None
		self.s2_pc_bit_count = None

	def getValue(self, yv):
		return getattr(self, yv)

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
		self.preprocess_time = None
		self.su_transmit_power = None
		self.rand_seed = None
		self.time_per_request = None
		self.secure_write_time = None
		self.path_loss = None

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
			elif isinstance(getattr(self, k), list):
				if len(getattr(self,k)) == 0:
					x = (k, None, 'list(float)')
				elif isinstance(getattr(self, k)[0], int):
					x = (k, None, 'list(int)')
				elif isinstance(getattr(self, k)[0], float):
					x = (k, None, 'list(float)')
				elif isinstance(getattr(self, k)[0], str):
					x = (k, None, 'list(str)')
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
		if var_name not in vars(self) or getattr(self, var_name) == None:
				return 'None'

		if var_order == None:
			if isinstance(getattr(self, var_name), list):
				if len(getattr(self, var_name)) > 0:
					return ','.join(map(str, getattr(self, var_name)))
				else:
					return ','
			else:
				return str(getattr(self, var_name))
		
		return ','.join(':'.join(map(str, vs)) for vs in zip(*tuple(getattr(self, var_name)[label] for label in var_order)))

	def hasAttr(self, yv):
		if isinstance(yv, str):
			return yv in vars(self)
		if isinstance(yv, tuple) and len(yv) == 2:
			return yv[0] in vars(self) and yv[1] in getattr(self, yv[0])

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

def runExperiment(param, no_run = False):
	args = ['../s2pc_stable', '-brief_out']

	var_flag_names = {
			'rand_seed': 'rand_seed', 'skip_s2pc': 'skip_s2pc',
			'central_entities': 'ces',
			'num_pu': 'npu', 'num_ss': 'nss', 'num_su': 'nsu', 'location_range': 'lr',
			'num_pr_per_pu': 'npr', 'pr_range': 'prr', 'out_filename': 'out',
			'unit_type': 'ut', 'propagation_model': 'pm',
			'ld_path_loss0': 'ld_pl0', 'ld_dist0': 'ld_d0', 'ld_gamma': 'ld_g',
			'splat_cmd': 'splat_cmd', 'ref_lat': 'ref_lat', 'ref_long': 'ref_long',
			'splat_dir': 'splat_dir', 'sdf_dir': 'sdf_dir', 'return_dir': 'return_dir',
			'num_ss_selection': 'nss_s', 'num_pu_selection': 'npu_s',
			'do_plaintext_split': 'do_pt_split', 'no_pr_thresh_update': '-no_pr_up',
			'selection_algo': 'sel_algo', 'secure_write_algo': 'sec_write_algo',
			'grid_x': 'grid_x', 'grid_y': 'grid_y',
			'ss_receive_power_alpha': 'rpa', 'ss_path_loss_alpha': 'pla',
			'num_float_bits': 'float_bits', 's2_pc_bit_count': 'bit_count'}

	# Check that all attrs have flags
	vs = []
	for attr in vars(param):
		if attr not in var_flag_names:
			raise Exception('No flag for attr ' + attr)

		vs.append((getattr(param, attr), var_flag_names[attr]))

	for val, flag in vs:
		if val != None:
			if isinstance(val, bool):
				if val == True:
					args += ['-' + flag]
			else:
				args += ['-' + flag, str(val)]

	if no_run:
		print ' '.join(args)
	else:
		rv = None
		for i in range(5):
			if args[0] != '../s2pc':
				print bcolors.WARNING + bcolors.UNDERLINE + 'Not running the most up to date version of code' + bcolors.ENDC

			process = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
			process.wait()

			if process.returncode != 0:
				if process.returncode == -9:
					print bcolors.FAIL + bcolors.BOLD + bcolors.UNDERLINE + 'Experiment Killed' + bcolors.ENDC
				else:
					print bcolors.FAIL + bcolors.BOLD + bcolors.UNDERLINE + 'Experiment failed with returncode: ' + str(process.returncode) + bcolors.ENDC

			# If there is any output to stderr, its assumed that something went wrong, and the experiment must be run again.
			err_out = [v for v in process.stderr]
			if len(err_out) != 0:
				print 'Encountered error:\n', ''.join(err_out)

			if process.returncode != 0 or len(err_out) != 0:
				continue

			rv = makeExperimentResult([v[:-1] for v in process.stdout])
			break

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
	elif isinstance(value, bool):
		type_str = 'bool'
		value_str = str(value)
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
		if isinstance(value[0], bool):
			type_str = 'list(bool)'
			value_str = ','.join(map(str, value))
		elif isinstance(value[0], int):
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
	if isinstance(value, bool):
		return 'bool'
	if isinstance(value, int):
		return 'int'
	if isinstance(value, float):
		return 'float'
	if isinstance(value, str):
		return 'str'
	raise Exception('Can\'t determine the type of: ' + str(value))

def convert(var_name, val_type, val_str):
	if val_type == 'None' or val_str == 'None':
		return var_name, None
	if val_type == 'bool':
		return var_name, val_str == 'True'
	if val_type == 'int':
		return var_name, int(val_str)
	if val_type == 'float':
		return var_name, float(val_str)
	if val_type == 'str':
		return var_name, val_str
	if val_type == 'list(bool)':
		return var_name, list(map(bool, val_str.split(',')))
	if val_type == 'list(int)':
		return var_name, list(map(int, val_str.split(',')))
	if val_type == 'list(float)':
		if val_str == '' or val_str == ',':
			return var_name, []
		return var_name, list(map(float, val_str.split(',')))
	if val_type == 'list(str)':
		return var_name, list(map(str, val_str.split(',')))
	if val_type == 'list(float,float,float)' or val_type == 'list(float,float)':
		simple_var_name = var_name[:var_name.find('(')]
		val_labels = var_name[var_name.find('(') + 1 : var_name.find(')')].split(',')

		vals = map(lambda x: tuple(map(float, x.split(':'))), val_str.split(','))
		vals = {label:[v[i] for v in vals] for i, label in enumerate(val_labels)}
		return simple_var_name, vals

	raise Exception('Unknown type: (' + str(var_name) + ', ' + str(val_type) + ', ' + str(val_str) + ')')

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Runs experiments. Has the parameters for each experiment saved, so they don\'t have to be reentered')
	
	current_time = datetime.now()

	parser.add_argument('-exp', '--experiment_identifier', metavar = 'EXPERIMENT_ID', type = str, nargs = '+', default = [], help = 'List of experiments to run. Must be value in (' + ', '.join(ALL_EXPERIMENT_IDS) + ')')
	parser.add_argument('-name', '--name_of_experiment', metavar = 'NAME', type = str, nargs = 1, default = [None], help = 'Includes the name in the output file')
	parser.add_argument('-out', '--out_file', metavar = 'OUT_FILE', type = str, nargs = 1, default = [None], help = 'File to write data to')
	parser.add_argument('-nt', '--num_tests', metavar = 'NUM_TESTS', type = int, nargs = 1, default = [1], help = 'Number of tests to run')
	parser.add_argument('-skip', '--skip_s2pc', action = 'store_true', help = 'If given, skips the s2pc algorithm')
	parser.add_argument('-hd', '--use_splat_hd', action = 'store_true', help = 'If given, uses splat-hd instead of splat')

	parser.add_argument('-nr', '--no_run', action = 'store_true', help = 'If given, doesn\'t run the tests, simply outputs the commands')

	args = parser.parse_args()
	experiments = args.experiment_identifier

	changes = []

	for experiment in experiments:
		if experiment == TEST:
			changes.append({NUM_SS_SELECTION: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 'num_pu_selection': [2], ('grid_x', 'grid_y', 'selection_algo'): [(100, 100, 'none')],
							'propagation_model': ['log_distance'], 'ld_path_loss0': [50], 'ld_dist0': [20], 'ld_gamma': [0.5],
							'num_pu': [10], PL_ALPHA: [2], RP_ALPHA: [2],
							'location_range': [100.0], 'num_ss': [100], 'num_su': [1], 'unit_type': ['db'],
							'num_pr_per_pu': [1], 'pr_range': [0.0]})
		if experiment == VARY_NUM_SS_SELECT:
			changes.append({NUM_SS_SELECTION: [1, 25, 50, 75, 100], ('grid_x', 'grid_y', 'selection_algo'): [(100, 100, 'none'), (250, 250, 'none'), (500, 500, 'none')],
							'propagation_model': ['longley_rice'],
							'num_pu': [10], PL_ALPHA: [2], RP_ALPHA: [2], 'location_range': [1000.0], 'num_ss': [500], 'unit_type': ['db'],
							'num_su': [100]})
		if experiment == TIMING_TEST:
			num_ss_s_test = {NUM_SS_SELECTION: [1, 10, 25, 50], 'num_pu_selection': [25], ('grid_x', 'grid_y'): [(1000, 1000)],
					'num_pr_per_pu' : [10], 'pr_range': [100.0],
					'propagation_model': ['log_distance'], 'ld_path_loss0': [50], 'ld_dist0': [20], 'ld_gamma': [0.5],
					'num_pu': [400], 'num_ss': [4000], 'num_su': [2],
					PL_ALPHA: [2], RP_ALPHA: [2], 'location_range': [10.0 * 1000.0], 'unit_type': ['db'],
					'secure_write_algo':['proposed'], 'central_entities': ['two_sms', 'sm_ks']}
			
			num_pu_s_test = {NUM_SS_SELECTION: [25], 'num_pu_selection': [1, 10, 25, 50], ('grid_x', 'grid_y'): [(1000, 1000)],
					'num_pr_per_pu' : [10], 'pr_range': [100.0],
					'propagation_model': ['log_distance'], 'ld_path_loss0': [50], 'ld_dist0': [20], 'ld_gamma': [0.5],
					'num_pu': [400], 'num_ss': [4000], 'num_su': [2],
					PL_ALPHA: [2], RP_ALPHA: [2], 'location_range': [10.0 * 1000.0], 'unit_type': ['db'],
					'secure_write_algo':['proposed'], 'central_entities': ['two_sms', 'sm_ks']}
			
			num_bits_test = {'s2_pc_bit_count': [64, 48, 32], NUM_SS_SELECTION: [25], 'num_pu_selection': [25], ('grid_x', 'grid_y'): [(1000, 1000)],
					'num_pr_per_pu' : [10], 'pr_range': [100.0],
					'propagation_model': ['log_distance'], 'ld_path_loss0': [50], 'ld_dist0': [20], 'ld_gamma': [0.5],
					'num_pu': [400], 'num_ss': [4000], 'num_su': [2],
					PL_ALPHA: [2], RP_ALPHA: [2], 'location_range': [10.0 * 1000.0], 'unit_type': ['db'],
					'secure_write_algo':['proposed'], 'central_entities': ['two_sms', 'sm_ks']}

			# secure_read_algo_test

			secure_write_algo_test = {NUM_SS_SELECTION: [25], 'num_pu_selection': [1, 10, 25, 50], ('grid_x', 'grid_y'): [(1000, 1000)],
					'num_pr_per_pu' : [10], 'pr_range': [100.0],
					'propagation_model': ['log_distance'], 'ld_path_loss0': [50], 'ld_dist0': [20], 'ld_gamma': [0.5],
					'num_pu': [400], 'num_ss': [4000], 'num_su': [2],
					PL_ALPHA: [2], RP_ALPHA: [2], 'location_range': [10.0 * 1000.0], 'unit_type': ['db'],
					'secure_write_algo':['proposed', 'spc'], 'central_entities': ['two_sms', 'sm_ks']}

			changes += [num_ss_s_test, num_pu_s_test, num_bits_test, secure_write_algo_test]

			for i in range(len(changes)):
				changes[i]['central_entities'] = ['two_sms']

		if experiment == PATH_LOSS_TEST:
			changes.append({NUM_SS_SELECTION: [1, 10, 25, 50], 'num_pu_selection': [25], ('grid_x', 'grid_y'): [(1000, 1000)],
							'propagation_model': ['log_distance'], 'ld_path_loss0': [50], 'ld_dist0': [20], 'ld_gamma': [0.5],
							'num_pu': [400], 'num_ss': [4000], 'num_su': [100],
							PL_ALPHA: [1, 2], RP_ALPHA: [1, 2], 'location_range': [10.0 * 1000.0], 'unit_type': ['db', 'abs'],
							'skip_s2pc': [True]})
		if experiment == MED_TEST:
			changes.append({NUM_SS_SELECTION: [1, 25, 50, 75, 100], 'num_pu_selection': [1, 5, 10],
							('num_pu', 'num_ss', 'location_range', 'grid_x', 'grid_y'):
								[(int(math.ceil(400 / (10000.0 * 10000.0) * x * x)),
								int(math.ceil(4000 / (10000.0 * 10000.0) * x * x)),
								x,
								int(x / 10),
								int(x / 10)) for x in (2000.0,)],
							'propagation_model': ['longley_rice'],
							PL_ALPHA: [2], RP_ALPHA: [2], 'unit_type': ['db']})
		if experiment == LARGE_TEST:
			changes.append({NUM_SS_SELECTION: [1, 10, 25, 50], 'num_pu_selection': [1, 10, 25,  50], ('grid_x', 'grid_y'): [(1000, 1000)],
							'num_pr_per_pu' : [10], 'pr_range': [100.0],
							'propagation_model': ['log_distance'], 'ld_path_loss0': [50], 'ld_dist0': [20], 'ld_gamma': [0.5],
							'num_pu': [400], 'num_ss': [4000], 'num_su': [10],
							PL_ALPHA: [2], RP_ALPHA: [2], 'location_range': [10.0 * 1000.0], 'unit_type': ['db']})
		if experiment == TEST_OUTPUT:
			changes.append({'num_pu': [10], 'num_ss': [100], 'num_su': [10], 'num_pr_per_pu' : [2], 'pr_range': [10.0],
							'location_range': [100.0], 'unit_type': ['db'],
							'propagation_model': ['single_lr'],
							'out_filename': ['../gen_out/test.txt']})
		if experiment == OUTPUT_RUN:
			changes.append({'num_pu': [400], 'num_ss': [4000], 'num_su': [1000], 'num_pr_per_pu' : [10], 'pr_range': [100.0],
							'location_range': [10.0 * 1000.0], 'unit_type': ['db'],
							'propagation_model': ['single_lr'],
							'out_filename': ['../gen_out/data1.txt']})

	num_experiments = sum(reduce(mul, [len(vals) for _, vals in change.iteritems()]) for change in changes) * args.num_tests[0]

	if num_experiments <= 0:
		raise Exception('No tests specified')

	out_values = []
	exp_num = 1
	durs = []
	begin_time = None
	for i in range(args.num_tests[0]):
		for change in changes:
			ns = change.keys()
			for vs in itertools.product(*tuple(vs for k, vs in change.iteritems())):
				print 'Running experiment', exp_num, 'of', num_experiments, 'with parameters:', ns, vs
				
				param = makeDefaultExperimentParam()

				if args.skip_s2pc:
					param.skip_s2pc = True

				if args.use_splat_hd:
					param.splat_cmd = 'splat-hd'

				for n, v in zip(ns, vs):
					if isinstance(n, tuple):
						for a, b in zip(n, v):
							vars(param)[a] = b
					else:
						vars(param)[n] = v

				start = datetime.now()
				if begin_time == None:
					begin_time = start
				result = runExperiment(param, no_run = args.no_run)
				if result != None:
					out_values.append((param, result))
				else:
					print bcolors.BOLD + bcolors.FAIL + 'Experiment failed to finish' + bcolors.ENDC
				end = datetime.now()
				durs.append(end - start)

				print 'Experiment took', bcolors.BOLD + bcolors.OKBLUE + str(end - start) + bcolors.ENDC

				if num_experiments - exp_num > 0:
					avg_dur = sum(durs, timedelta(0)) / len(durs)
					expected_end = datetime.now() + avg_dur * (num_experiments - exp_num)
					if expected_end.date() == datetime.now().date():
						expected_end_str = expected_end.strftime('%-I:%M:%S.%f %p')
					else:
						expected_end_str = expected_end.strftime('%A %b %-d %-I:%M:%S.%f %p')

					print 'Expected end at', bcolors.BOLD + bcolors.OKGREEN + expected_end_str + bcolors.ENDC
				exp_num += 1

	end = datetime.now()
	print 'All Experiments took', bcolors.BOLD + bcolors.WARNING + str(end - begin_time) + bcolors.ENDC

	if len(out_values) > 0:
		out_file = args.out_file[0]
		if out_file == None:
			out_file = current_time.strftime('data/out_%b-%d_%H:%M:%S_' + (args.name_of_experiment[0] if args.name_of_experiment[0] != None else '-'.join(experiments)) + '.txt')

		writeOutParamResult(out_file, out_values)
	
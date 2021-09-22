import os.path 
import numpy as np 
from sklearn.neighbors import KernelDensity as KD
from scipy.interpolate import CubicSpline
from scipy import integrate
from scipy.optimize import fmin, root
import matplotlib.pyplot as plt
from statistics import median

def parse_config_file(config_path = os.path.abspath('.'), config_filename = 'analysis_config.txt'):
    '''
    Reads in the configuration file 'analysis_config.txt' and converts the specified data to the appropriate type. 
    '''
    config_dict = default_config()
    config_file_path = os.path.join(config_path, config_filename)

    with open(config_file_path, 'r') as file:
        line_count = 0
        for line in file:
            line_count += 1
            if line[0]=="#" or not line.strip(): continue
            (key, val) = line.split('=')
            key = key.strip()
            val = val.strip().strip('"').strip("'")
            if key in ['PLOT_ERROR_BARS', 'JITTER', 'PLOT_LINE', 'PLOT_DATA_POINTS', 'PLOT_CURVE_ERROR', 'REL_POT_TO_REF']: 
                val = True if "true" in val.lower() else False
            elif key in ['LC_VALUES']: 
                val = np.array([float(i.strip()) for i in val.split(",")])/100.
            elif key in ['RHO', 'BETA_PRIOR', 'BETA_PRIOR_0', 'JITTER_FACTOR', 'ALPHA', 'CURVE_CI', 'LC_CI', 'ERROR_BAR_CI', 'REL_POT_CI']: 
                val = float(val)/100. if key in ['CURVE_CI', 'LC_CI', 'REL_POT_CI', 'ERROR_BAR_CI'] else float(val)
                if key == 'JITTER_FACTOR':
                    if val < 0: val = 0
                    elif val > 0.5: val = 0.4
            elif key in ['BOOTSTRAP_ITERS', 'N_POINTS', 'NCPU', 'ERROR_BAR_CUTOFF']: 
                val = int(val)
            config_dict[key] = val
    return config_dict

def default_config():
    config_dict = {
                    'ARCHIVE_DIR': None,
                    'OUT_DIR': None,
                    'NCPU': -1,
                    'BOOTSTRAP_ITERS': 1000,
                    'RHO': 0.10,
                    'BETA_PRIOR': 0.0,
                    'BETA_PRIOR_0': 0.25,
                    'CI_METHOD': "HPDI",
                    'CURVE_TYPE': "best",
                    'FIT_METHOD': "BFGS",
                    'N_POINTS': 151,
                    'PLOT_LINE': True,
                    'LINE_COLOR': 'black',
                    'LINE_STYLE': 'solid',
                    'PLOT_DATA_POINTS': True,
                    'POINT_COLOR': 'black',
                    'MARKER_STYLE': '.',
                    'PLOT_ERROR_BARS': True,
                    'ERROR_BAR_CI': 80,
                    'JITTER': True,
                    'JITTER_FACTOR': 0.1,
                    'PLOT_CURVE_ERROR': True,
                    'ALPHA': 0.5,
                    'CURVE_CI_COLOR': 'lightgrey',
                    'CURVE_CI': 95 ,
                    'LC_VALUES': [50, 90],
                    'LC_CI': 95,
                    'REL_POT_CI': 95,
                    'REFERENCE_COMPOUND': 'malathion',
                    'REL_POT_TO_REF': True, 
                    'ERROR_BAR_CUTOFF': 10
                    }
    return config_dict

def fix_num_output(func):
    def fix_scinot(*args, **kwargs):
        return func(*args,**kwargs).replace("e+0","e").replace("E+0","e").replace("e-0","e-").replace("E-0","e-").replace("e+","e")
    return fix_scinot

def check_library_change(cmpd_options, dict_options):
    changed = False
    check_list = ["BOOTSTRAP_ITERS", "RHO", "BETA_PRIOR", "CURVE_TYPE"]
    for item in check_list:
        if cmpd_options[item] != dict_options[item]:
            changed = True
    return changed

def calc_HPDI_CI(kernel_density = None, kernel_sample = None, CI_level = 0.95, n_samples = 500):
    alpha = 1 - CI_level
    lb = alpha/2
    ub = 1- lb
    if kernel_density is None and kernel_sample is not None:
        kernel_density = KD(kernel='gaussian', bandwidth=0.5).fit(kernel_sample)
    # if kernel_sample 
    if kernel_density is not None:
        x = kernel_density.sample(100000)
        med = median(x)
        minx, maxx = np.quantile(x, [0.01, 0.99], interpolation='linear', axis = 0)
        print(lb, ub)
        l1, l2, l3 = np.quantile(x, [lb, 0.5, ub], interpolation='linear', axis = 0)
        ranx = max(x) - min(x)
        if ranx > 15: #Too big to give decent response curve.
            return np.quantile(x, [lb, 0.5, ub], interpolation='linear', axis = 0)
        else:
            lb, ub = minx - 1, maxx + 1
            x =  np.linspace(lb, ub, round(float(ranx)*15+1))
            y = kernel_density.score_samples(x)
            y = np.exp(y)
            spline = CubicSpline(x.squeeze(), y)
            

            # xm = kernel_density.sample(200)
            # fig, ax = plt.subplots()
            # plt.plot(xm,  np.exp(kernel_density.score_samples(xm)), '.g')
            # plt.plot(xm, spline(xm), '.r')
            # plt.show()

    print("finding min")

    def errfn(p, CI_level, spline, lb, ub):
        # print(p)
        def fn( x ): 
            # print(x, spline(x), p, spline(x) if spline(x) > p else 0)
            return spline(x) if spline(x) > p else 0

        prob = integrate.quad(fn, lb, ub,limit = 100)[0]
        return (prob - CI_level)**2

    p = fmin(errfn, x0=0.1, args=(CI_level, spline, lb, ub))[0]
    
    def root_finders(x, spline, p): return spline(x) - p 
    minval = root(root_finders, x0=lb, args=(spline, p)).x
    maxval = root(root_finders, x0=ub, args=(spline, p)).x
    print(np.array([l1, l2, l3]) )
    print(np.array([minval, med, maxval]) )
    return np.array([minval, med, maxval]) 

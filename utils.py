import os.path 
import numpy as np 
from sklearn.neighbors import KernelDensity as KD
from scipy.interpolate import CubicSpline
from scipy import integrate
from scipy.optimize import fmin, root
import matplotlib.pyplot as plt
from statistics import median

ET_VARS = ["equal-tailed", "equal", "et", "even", "even-tailed"]

def parse_config_file(config_path = os.path.abspath('.'), config_filename = 'analysis_config.txt', **kwargs):
    '''
    Reads in the configuration file 'analysis_config.txt' and converts the specified data to the appropriate type. 
    '''
    config_dict = default_config()
    config_file_path = os.path.join(config_path, config_filename)

    #first read config file
    with open(config_file_path, 'r') as file:
        line_count = 0
        for line in file:
            line_count += 1
            if line[0]=="#" or not line.strip(): continue
            line = line.split('#')[0] #don't read anything after the first '#'
            (key, val) = line.split('=')
            key, val = config_parse_helper(key, val)
            #default behaviour: val is simply the stripped string. 
            config_dict[key] = val

    #then use kwargs, if present. NB: This supersedes the config file. 
    for key, val in **kwargs.items():
        key, val = config_parse_helper(key, val)
        config_dict[key] = val

    return config_dict


def config_parse_helper(key, val):
    key = key.strip()
    val = val.strip().strip('"').strip("'")
    if key in ['PLOT_ERROR_BARS', 'JITTER', 'PLOT_LINE', 
                'PLOT_DATA_POINTS', 'PLOT_CURVE_ERROR', 'REL_POT_TO_REF']: 
        val = True if "true" in val.lower() else False
    elif key in ['LC_VALUES']: 
        val = np.array([float(i.strip()) for i in val.split(",")])/100.
    elif key in ['RHO', 'BETA_PRIOR', 'BETA_PRIOR_0', 
                'JITTER_FACTOR', 'ALPHA', 'CURVE_CI', 
                'LC_CI', 'ERROR_BAR_CI', 'REL_POT_CI', 
                'EXTRAPOLATION_FACTOR']: 
        val = float(val)/100. if key in ['CURVE_CI', 'LC_CI', 'REL_POT_CI', 'ERROR_BAR_CI'] else float(val)
        if key == 'JITTER_FACTOR':
            if val < 0: val = 0
            elif val > 0.5: val = 0.4
    elif key in ['BOOTSTRAP_ITERS', 'N_POINTS', 'NCPU', 'ERROR_BAR_CUTOFF']: 
        val = int(val)
    elif key in ["CI_METHOD"]:
        val = val.lower()
    elif key in ["OUTPUT_PDF_NAME"]:
        if val[-4:].lower() in ['.pdf', '.tex']: val = val[:-4] #strip pdf ending

def default_config():
    config_dict = {
                    'ARCHIVE_PATH': None,
                    'SAVE_PATH': None,
                    'OUTPUT_CSV_NAME': 'output.csv',
                    'OUTPUT_PDF_NAME': 'graphs.pdf',
                    'INPUT_PATH': None,
                    'MASK_RCNN_SAVE': True,
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
                    'ERROR_BAR_CUTOFF': 10,
                    "EXTRAPOLATION_FACTOR": 2
                    }
    return config_dict

def fix_num_output(func, *args, **kwargs):
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

# def calc_multimodal_HPDI_CI(kernel_density = None, kernel_sample = None, CI_level = 0.95, n_samples = 500):
#     '''
#     Determines HDPI CIs for multimodal distributions. Briefly, the kernel density function K(x)
#     is approximated by a cubic spline, and a y-value is found at which the integral of the function
#     [K`(x) = K(x) if K(x) > y else 0] is equal to 1-alpha. 
#     '''
#     alpha = 1 - CI_level
#     lb = alpha/2
#     ub = 1- lb
#     if kernel_density is None and kernel_sample is not None:
#         kernel_density = KD(kernel='gaussian', bandwidth=0.5).fit(kernel_sample)
#     # if kernel_sample 
#     if kernel_density is not None:
#         x = kernel_density.sample(100000)
#         med = median(x)
#         minx, maxx = np.quantile(x, [0.01, 0.99], interpolation='linear', axis = 0)
#         l1, l2, l3 = np.quantile(x, [lb, 0.5, ub], interpolation='linear', axis = 0)
#         ranx = max(x) - min(x)
#         if ranx > 15: #Too big to give decent response curve.
#             return np.quantile(x, [lb, 0.5, ub], interpolation='linear', axis = 0)
#         else:
#             lb, ub = minx - 1, maxx + 1
#             x =  np.linspace(lb, ub, round(float(ranx)*15+1))
#             y = kernel_density.score_samples(x)
#             y = np.exp(y)
#             spline = CubicSpline(x.squeeze(), y)

#     def errfn(p, CI_level, spline, lb, ub):
#         # print(p)
#         def fn( x ): return spline(x) if spline(x) > p else 0

#         prob = integrate.quad(fn, lb, ub,limit = 100)[0]
#         return (prob - CI_level)**2

#     p = fmin(errfn, x0=0.1, args=(CI_level, spline, lb, ub))[0]
    
#     def root_finders(x, spline, p): return spline(x) - p 
#     minval = root(root_finders, x0=lb, args=(spline, p)).x
#     maxval = root(root_finders, x0=ub, args=(spline, p)).x
#     return np.array([minval, med, maxval]) 

def CI_helper(kernel_sample, CI_level = 0.95, min_sample_size = 100000, resample = True):
    
    def est_gaussian_kernel_bandwidth(data):
        ''' 
        Uses Scott’s rule of thumb with IQR adjustment
        According to the wikipedia page https://en.wikipedia.org/wiki/Kernel_density_estimation,
        which cites
        Silverman, B.W. (1986). Density Estimation for Statistics and Data Analysis. 
            London: Chapman & Hall/CRC. p. 45. ISBN 978-0-412-24620-3.
        '''
        IQR = np.quantile(kernel_sample, [0.25, 0.75], interpolation='linear', axis = 0)
        mid = (IQR[1] - IQR[0])/1.34
        std = np.std(data, axis = 0)
        mult = np.minimum(mid, std)
        bw = 0.9*mult*(len(data)**-0.2)
        return bw

    n = len(kernel_sample)
    if n < min_sample_size and resample:
        if kernel_sample.ndim == 1 : kernel_sample = kernel_sample.reshape(-1, 1) #reshape from weird shape
        bw = est_gaussian_kernel_bandwidth(kernel_sample)
        for i in range(len(bw)):
            #again reshape from weird shape
            kernel_density = KD(kernel='gaussian', bandwidth=bw[i]).fit(kernel_sample[:,i].reshape(-1, 1))
            if i == 0: sample = kernel_density.sample(min_sample_size)
            else: sample = np.append(sample, kernel_density.sample(min_sample_size), axis = 1)
    else: sample = kernel_sample.squeeze()
    
    return np.sort(sample, axis= 0)

def calc_HPDI_CI(kernel_sample, CI_level = 0.95, min_sample_size = 100000, resample = True):
    """
    Internal method to determine the minimum interval of a given width
    """
    kernel_sample = CI_helper(kernel_sample, 
        CI_level = CI_level, 
        min_sample_size = min_sample_size, 
        resample = resample)

    n = len(kernel_sample)
    med = np.median(kernel_sample, axis = 0)
    # print(med)

    interval_idx_inc = int(np.floor(CI_level*n))
    n_intervals = n - interval_idx_inc

    interval_width = kernel_sample[interval_idx_inc:] - kernel_sample[:n_intervals]
    if len(interval_width) == 0: raise ValueError('Too few elements for interval calculation')

    min_idx = np.argmin(interval_width, axis=0)
    hdi_min = []
    hdi_max = []
    col = 0
    if interval_width.ndim > 1:
        for idx in min_idx:
            # print(col, idx)
            hdi_min.append(kernel_sample[idx,col])
            hdi_max.append(kernel_sample[idx+interval_idx_inc,col])
            col += 1
    else:
        hdi_min = kernel_sample[min_idx]
        hdi_max = kernel_sample[min_idx+interval_idx_inc]
    return np.array([hdi_min, med, hdi_max]) 

def calc_ET_CI(kernel_sample, CI_level = 0.95, min_sample_size = 25000, resample = True):
    """
    Internal method to determine the minimum interval of a given width
    """
    kernel_sample = CI_helper(kernel_sample, 
        CI_level = CI_level, 
        min_sample_size = min_sample_size, 
        resample = resample)

    alpha = 1 - CI_level
    lb = alpha/2
    ub = 1- lb
    
    return np.quantile(kernel_sample, [lb, 0.5, ub], interpolation='linear', axis = 0)

@fix_num_output
def format_lc_val(val):
    # print(val)
    if val < 1e-2 or val >1e3: 
        return f"{val:.2e}"
    else: return f"{val:.3g}"

def CI_to_string(ci1, ci2):
    return '"[' + format_lc_val(ci1) +", "+ format_lc_val(ci2) + ']"'

def format_LC_to_CSV(CI_array):
    def append_line(output, line):
        output.append(format_lc_val(line[1]))
        output.append(CI_to_string(line[0],line[2]))
    output = []
    if CI_array.ndim > 1:
        for line in CI_array: append_line(output, line)
    else: append_line(output, CI_array)
    return output


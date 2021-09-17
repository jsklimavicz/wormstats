import os.path 
import numpy as np 


def parse_config_file(config_path = os.path.abspath('.'), config_filename = 'analysis_config.txt'):
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
            if key in ['ERROR_BARS', 'JITTER', 'PLOT_LINE', 'PLOT_DATA_POINTS', 'PLOT_CURVE_ERROR']: 
                val = True if "true" in val.lower() else False
            elif key in ['LC_VALUES']: 
                val = np.array([float(i.strip()) for i in val.split(",")])/100.
            elif key in ['RHO', 'SCALE', 'JITTER_FACTOR', 'ALPHA', 'CURVE_CI', 'LC_CI', 'ERROR_BAR_CI']: 
                val = float(val)/100. if key in ['CURVE_CI', 'LC_CI', 'ERROR_BAR_CI'] else float(val)
            elif key in ['BOOTSTRAP_ITERS', 'N_POINTS']: 
                val = int(val)
            config_dict[key] = val
    return config_dict

def default_config():
    config_dict = {
                    'ARCHIVE_DIR': None,
                    'OUT_DIR': None,
                    'BOOTSTRAP_ITERS': 1000,
                    'RHO': 0.10,
                    'SCALE': 0.10,
                    'CI_METHOD': "HPDR",
                    'CURVE_TYPE': "best",
                    'FIT_METHOD': "BFGS",
                    'N_POINTS': 151,
                    'PLOT_LINE': True,
                    'LINE_COLOR': 'black',
                    'LINE_STYLE': 'solid',
                    'PLOT_DATA_POINTS': True,
                    'POINT_COLOR': 'black',
                    'MARKER_STYLE': '.',
                    'ERROR_BARS': True,
                    'ERROR_BAR_CI': 80,
                    'JITTER': True,
                    'JITTER_FACTOR': 0.1,
                    'PLOT_CURVE_ERROR': True,
                    'ALPHA': 0.5,
                    'CURVE_CI_COLOR': 'lightgrey',
                    'CURVE_CI': 95 ,
                    'LC_VALUES': [50, 90],
                    'LC_CI': 95,
                    'REFERENCE_COMPOUND': 'malathion',
                    }
    return config_dict


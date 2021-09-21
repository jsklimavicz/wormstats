# wormstats

## Intro

This set of files is to calculate Bayesian log-likelihood dose-response curves based on bioassay live-dead data. Briefly, live/dead count data, associated with compounds/codes, dates, and plate rows and columns, is read into the program. These live/dead counts are then used to create bootstrapped samples of live/dead ratios using a beta distribution, and dose-response curves are fit to this bootstrapped data. These curves are then used to generate credible intervals for the best fit curve, which is taken as the median response at each concentration for which the curve is calculated. 

The output data is in .csv value, and includes estimates and credible intervals of user-requested LC values. The user can also specify a reference/control compound, and estimates and credible intervals of potency ratios of each compound to this reference compound are given. 

Images of each dose-response curve, the curve credible interval, the actual live/dead data, and the bayesian credible intervals for the actual data are produced. 

## Advanced Description



Because the numerical optimization of many dose-response curves is somewhat computationally expensive, the program saves fitted curve parameters by pickling. 


## User-Specified Values

The program accepts user-specified values in the file ``analysis_config.txt``. The allowed variables are broken down into groups roughly as follows.

### Saving Options

These options affect where/how/whether data is saved. Specifically, ``IMAGE_DIR`` specifies the directory in which images should be saved; the program will attempt to create the directory if it does not exist. The ``OUTPUT_DIR`` variable specifies the directory in which the output .csv file and composite graphs are saved. The ``ARCHIVE_DIR`` specifies where pickled data should be saved (see the Advanced Description), while ``INPUT_PATH`` is the full path to the .csv file with live/dead data.

### Graphing Options

#### Included Components 
The program allows the user to specify many parts of the graphing process. First, the user can specify whether certain parts of the graph are even included. The variables `PLOT_LINE`, `PLOT_DATA_POINTS`, `PLOT_CURVE_ERROR`, and `PLOT_ERROR_BARS` plot the best-fit line, the actual live/dead ratio data points, the credible interval of the curve, and the credible intervals for each data point, respectively. Each of these values can either be set to `True` (default for all) if they are to be included, or `False`, if they should not be plotted. 

#### Colors and Styles
The best-fit line, by default, is a black, solid line; however, the color can be changed by adjusting the value of `LINE_COLOR` to a [permitted value][1], while the line style can be adjusted by using a [permitted value][2] for the `LINE_STYLE`. Likewise, the style of the [marker][3] for each data point can be adusted through `MARKER_STYLE` (default: `.`, which specifies a point), and the [color][1] through `POINT_COLOR` (default: black). The credible interval lines for each data point are automatically set to the `POINT_COLOR`. 

The [color][1] of the credible region for the curve may be adjusted through `CURVE_CI_COLOR` (default: lightgray), while the transparency of the curve may be set through the variable `ALPHA` (default 0.5), which takes on values between 1 (fully opaque) and 0 (fully transparent). 

#### The Look

To plot the smoothest curves, it is necessary to have sufficient x-values included when calculating the curve. The variable `N_POINTS` (default: 151) sepcifies the number of points to use. The trade-off is that while larger values give nicer curves, smaller values result in smaller vector-image files and faster computational times. 

The values of the credible intervals for both the error bars and the dose-response curve can be set; the former through `ERROR_BAR_CI` (default: 80), the latter, `CURVE_CI` (default: 95). Both of these variables take values between 0 and 100, which specify the desired credible interval (e.g. `ERROR_BAR_CI` = 80 specifies that each error bar should represent the 80% credible interval of each data point). 

To prevent overcrowding of data points and error bars, the option to `JITTER` the data points is given. Setting this variable to `True` (default) mimics the behavior of `geom_jitter` in `ggplot2` in R, and randomly adjusts the x-value of the data. The amount by which the data can be jittered is given by `JITTER_FACTOR` (default: 0.15), which should be between 0 (no jitter) and 0.5 (jittering by more than this amount would allow data differing by a factor of 2 to overlap). Attempting to set `JITTER_FACTOR` below 0 or above 0.5 will automatically result in `JITTER_FACTOR` being set to 0 or 0.4, respectively. Additionally, the graph may get very crowded when many replicates are performed; in this instance, it may be desirable to omit the plotting of error bars for these compounds only. This behaviour can be controlled by setting `ERROR_BAR_CUTOFF` to the number of reps at which the error bars are no longer plotted (default: 10).

### Data-Fitting Options

Because the curve-fitting process is the most computationally-expensive component of the analysis, this portion has been written to be performed in parallel. The number of CPUs that should be used for computation should be set with `NCPU`. The default behavior is to use all available CPUs, which is automatically determined if `NCPU` = -1 (default). Alternatively, this variable may be set to any positive integer up to the number of available CPUs. 

`BOOTSTRAP_ITERS` (default: 2500) speficies the number of bootstrapped curves to calculate. Higher numbers of bootstrap iterations will produce more reliable confidence intervals and curves, but correspondingly, this will take more computational time. For best results, this value should be no lower than 500. 

As discussed in the Advanced Description, the bootstrapped live/dead probabilities are created using beta distributions, since the live count out of the total count is binomial in nature. Because the posterior distribution is a beta distribution, the user is given some choice of prior through the variable `BETA_PRIOR`. Beta variables are then drawn from the distribution Beta(live count + `BETA_PRIOR`, dead count + `BETA_PRIOR`). Setting `BETA_PRIOR` to 0.5 is equivalent to Jeffreys' prior, while using 0.0 (default) is roughly equivelent to Heldane's prior. To avoid the [sunrise paradox][4], the `BETA_PRIOR_0` value can be set to provide a seperate beta prior when the live count or dead count is equal to 0 (default: 0.25). To prevent undefined values, `BETA_PRIOR_0` takes on a minimum value of sqrt(machine_epsilon).

The `RHO` parameter is used to set correlation between the beta variables contained within a single replicate. See Advanced Description for more details. This paramter should be set to between -0.25 and roughly 0.6 (default: 0.10). 

`CI_METHOD` specifies the method of confidence interval calculation. The default behavior (`HPDI`) is to calculate the highest posterior density interval (this is also the optimal CI in terms of decision theory). Setting this value to ``equal-tailed`` will produce credible intervals with equal probabilities of being above or below the median, and equal tail weights. This behavior is set for the credible intervals for both graphing and for the output .csv file. 


### Analysis Options

The user may specify one or more desired values for LC estimations as a comma-seperated list using the variable `LC_VALUES`. The default values (`50, 90`) results in the estimation of an LC50 and an LC90, and the corresponding credible intervals. The credible interval size is specified using `LC_CI` (default: 95). Relative potencies can also be calculated with respect to a reference compound (defined using `REFERENCE_COMPOUND`, default: malathion). The width of the credible interval for the potency ratio is specified using `REL_POT_CI`.
 
## Files
This packages contains the following files:
- `merlinanalyzer.py`, which contains the class object that drives the data analysis.
- `compound.py` defines the class that contains all the data and calculated parameters for a compound.
- `curvell.py` does most of the mathematical calculations for curve fitting, parameter estimation, etc.
- `corr_beta.py` generates the correlated random beta variables for bootstrapping purposes. 
- `utils.py` contains several useful functions, including the function that creates the options dictionary from `analysis_config.txt`.
- `analysis_config.txt` is used to set the user-chosen parameters. 







[1]: https://matplotlib.org/stable/gallery/color/named_colors.htm
[2]: https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
[3]: https://matplotlib.org/stable/api/markers_api.html#module-matplotlib.markers
[4]: https://en.wikipedia.org/wiki/Sunrise_problem

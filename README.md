# wormstats

## Intro

This set of files is to calculate Bayesian log-likelihood dose-response curves based on bioassay live-dead data. Briefly, live/dead count data, associated with compounds/codes, dates, and plate rows and columns, is read into the program. These live/dead counts are then used to create bootstrapped samples of live/dead ratios using a beta distribution, and dose-response curves are fit to this bootstrapped data. These curves are then used to generate credible intervals for the best fit curve, which is taken as the median response at each concentration for which the curve is calculated. 

The output data is in .csv value, and includes estimates and credible intervals of user-requested LC values. The user can also specify a reference/control compound, and estimates and credible intervals of potency ratios of each compound to this reference compound are given. 

Images of each dose-response curve, the curve credible interval, the actual live/dead data, and the bayesian credible intervals for the actual data are produced. 

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



## Advanced Description

Because the numerical optimiztion of many dose-response curves is somewhat computationally expensive, the program saves fitted curve parameters by pickling. 

### Files








[1]: https://matplotlib.org/stable/gallery/color/named_colors.htm
[2]: https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
[3]: https://matplotlib.org/stable/api/markers_api.html#module-matplotlib.markers

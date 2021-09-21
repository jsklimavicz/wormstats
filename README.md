# wormstats

## Intro

This set of files is to calculate Bayesian log-likelihood dose-response curves based on bioassay live-dead data. Briefly, live/dead count data, associated with compounds/codes, dates, and plate rows and columns, is read into the program. These live/dead counts are then used to create bootstrapped samples of live/dead ratios using a beta distribution, and dose-response curves are fit to this bootstrapped data. These curves are then used to generate credible intervals for the best fit curve, which is taken as the median response at each concentration for which the curve is calculated. 

The output data is in .csv value, and includes estimates and credible intervals of user-requested LC values. The user can also specify a reference/control compound, and estimates and credible intervals of potency ratios of each compound to this reference compound are given. 

Images of each dose-response curve, the curve credible interval, the actual live/dead data, and the bayesian credible intervals for the actual data are produced.

## User-Specified Values

The program accepts user-specified values in the file <analysis_config.txt>. 

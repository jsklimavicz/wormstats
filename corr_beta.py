import math
import scipy.stats as st
import numpy as np


def corr_matrix(plate_ids, rho = 0.10):
	'''
	Generates a len(plate_ids) x len(plate_ids) block matrix with the number of 
	blocks equal to the numebr of unique plates for a compound. 
	Rho is a correlation coefficient for a Gaussian MVN distribution for generating
	correlated beta variables by means of a copula. 
	'''
	mat = np.eye(len(plate_ids))
	for p in set(plate_ids):
		locations = [j for j, x in enumerate(plate_ids) if x == p]
		for i in locations:
			for j in locations:
				if i > j:
					val = rho/(2**(abs(i-j)-1))
					mat[i,j] = val
					mat[j,i] = val
	return mat

def corr_beta_vars(plate_ids, live_count, dead_count, size = 1, scale = 0.5, rho = 0.10):
	#TODO: Implement Heldane's prior. 
	if scale < np.sqrt(np.finfo(float).eps): scale =  np.sqrt(np.finfo(float).eps)
	#Generates the correlated beta variables
	#sanity checks
	assert len(plate_ids) == len(live_count)
	assert len(dead_count) == len(live_count)
	#get a correlation matrix
	corr_mat = corr_matrix(plate_ids, rho = rho)
	#generate multivariate normal random variables
	norm_var = np.random.default_rng().multivariate_normal(np.zeros_like(plate_ids), corr_mat, size = size)
	#find normal cdf values 
	norm_prob = st.norm.cdf(norm_var)
	#generate beta variables from cdf values
	scale_add = np.where(np.logical_or(dead_count == 0, live_count == 0), 0.25, scale)
	beta_var = st.beta.ppf(norm_prob, dead_count + scale_add, live_count + scale_add)
	return beta_var

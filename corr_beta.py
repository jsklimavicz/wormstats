import math
import scipy.stats as st
import numpy as np


def corr_matrix(plate_ids, rho = 0.25):
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

def corr_beta_vars(plate_ids, live_count, dead_count, size = 1, scale = 0.5, rho = 0.25):
	assert len(plate_ids) == len(live_count)
	assert len(dead_count) == len(live_count)
	corr_mat = corr_matrix(plate_ids, rho = rho)
	norm_var = np.random.default_rng().multivariate_normal(np.zeros_like(plate_ids), corr_mat, size = size)
	norm_prob = st.norm.cdf(norm_var)
	beta_var = st.beta.cdf(norm_prob, live_count + scale, dead_count + scale)
	return beta_var

# Deconvolution-in-well-test-analysis

* __CHARM.R__ optimization step
* __DEzs.R__ posterior shape retrieval step
* __00c_sample_p0q.R__ retrieve the marginalized parameters (optional)
* __functions.R__ auxiliary functions
* __1a - 05b__ vizualization scripts

The first MCMC variation (CHARM) is acting as an optimization step, in order to find the MAP. The second MCMC variation (DEzs) is using information from the previous step in order to retrieve the posterior distribution.

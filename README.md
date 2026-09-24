# SMCM-SL-MH-Algorithm
3-cloud SMCM learning algorithm using Metropolis-Hastings paired with a Synthetic Likelihood measure with Summary Statistics designed specifically for the SMCM. Provided training Data set was produced using a  simple Single Column model coupled to the SMCM.

This is code and data used to produce the results in the submitted paper 'BAYESIAN SYNTHETIC LIKELIHOOD ALGORITHMS TO1
CALIBRATE A STOCHASTIC MODEL FOR CONVECTION' by Carlos Sevilla, Boualem Khouider and Carsten Abraham. For questions on this code or data please contact csevilla@uvic.ca. For information on the motivation, purpose and results of this program please see the named paper.

BSL_MCMC.m - The Main file, houses all inference settings, MH proposal distribution variance, initial parameter values. Whitening of the Summary Statistics Covariance Matrix is Available but commented out. 

equilibriumdistributionSL.m - Computes the cloud transition Rates Given predictor variables

flow.m - computes the summary statistics

gammabb.m - Computes negative exponential of predictive variables

logmvnpdf.m - Computes the log likelihood

Sample_Path.m - Produces sample paths of SMCM given training data set of predictive variables

SMC_Only.m - Computes cloud transitions given transition rates for a given time step

Whitening.m - Alternate main file to compute whitening matrix


The Timeseries directory contains resultant MCMC sample chains as reported in the above paper. 


# Program R Statistical Population Reconstruction of Moose

## Overview
This repository contains R code for statistical population reconstruction of moose (Alces alces) in Northeastern Minnesota. The reconstruction combines several data sources: age-at-harvest data from both state and tribal harvest, aerial survey estimates, and data from radio-collared animals. The reconstruction then creates estimates for population dynamics over time 

## Table of Contents
- [Set Up Work Space](#set-up-work-space)
- [Import Raw Data](#import-raw-data)
- [Synthesize Data](#synthesize-data)
- [Define Objective Function](#define-objective-function)
- [Perform Numerical Optimization](#perform-numerical-optimization)
- [Extract Uncertainty Estimates](#extract-uncertainty-estimates)

## Set Up Work Space
The code begins by clearing the global environment and importing the necessary R packages, including BB, pso, readxl, ggplot2, numDeriv, and truncnorm. These packages are essential for data manipulation, numerical optimization, and plotting.

## Import Raw Data
The next section imports the raw data, including age-at-harvest matrices for state and tribal harvests, aerial survey data, and telemetry data. This was the most current data at the time, and would need to be updated as more data becomes available.

## Synthesize Data
This section extracts and synthesizes data into the formate that is necessary to do maximum likelihood estimation - calculating the means and standard deviations of the collected data.

## Define Objective Function
The core of the code defines the objective function using a multinomial likelihood formulation - this creates several likelihood functions from the data sources and combines these into a single measure of how well the input parameters of harvest rates, survival rates, and abundance estimates fit the collected data. These likelihood functions are calculated to be log-likelihoods, meaning that the sum is taken to return the final measure.

## Perform Numerical Optimization
Numerical optimization is done using particle swarm optimization is used to find the best-fit model, and the results are stored in a data frame, including the estimates for the harvest rates, survival rates and population abundance estimates, model settings, and the AIC value.

## Extract Uncertainty Estimates
The final section extracts uncertainty estimates for the best-fit model, including standard errors for parameters and stochastic abundance estimates which provide an understanding of how reliable the results are. This was then later used to compare multiple models with different parameters fixed across years.


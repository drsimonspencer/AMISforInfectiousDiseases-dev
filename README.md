# AMIS for Infectious Diseases

[![R-CMD-check](https://github.com/drsimonspencer/AMISforInfectiousDiseases-dev/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/drsimonspencer/AMISforInfectiousDiseases-dev/actions/workflows/R-CMD-check.yaml)

## Description

This package provides an implementation of the Adaptive Multiple
Importance Sampling algorithm, as described in 

Integrating geostatistical maps and infectious disease transmission models using adaptive multiple importance sampling.
Renata Retkute, Panayiota Touloupou, María-Gloria Basáñez, T. Deirdre Hollingsworth and Simon E.F. Spencer (2021).
_Annals of Applied Statistics_, 15 (4), 1980-1998. DOI https://doi.org/10.1214/21-AOAS1486 

This repository contains the development version of the package.

## Installation

Make sure you have the package [devtools](https://devtools.r-lib.org/)
installed. Then

```R
devtools::install_github("drsimonspencer/AMISforInfectiousDiseases-dev")
```

## Usage

`amis` is the main function of the package. It takes a
geostatistical map, a transmission model and a prior distribution for parameters, 
and returns a list of objects that includes sampled parameters and their associated weights 
in each location at each time point.

```R
amis_output <- AMISforInfectiousDiseases::amis(prevalence_map, transmission_model, prior, amis_params)
```

- `prevalence_map`: A matrix representing the geostatistical map, with one
  row per pixel (if there is one time point); or a list with matrices, 
  each one representing the geostatistical map for a time point.
- `transmission_model`: A function implementing the model. It can be anything, as
  long as it conforms to a specific interface. See defining a model
  function.
- `prior`: A list containing the functions dprior and rprior (density and RNG, respectively).
- `amis_params`: A list containing the parameters for the AMIS algorithm, such as:
  - `n_samples`: The number of sample parameters to draw at each iteration.
  - `target_ess`: The target effective sample size.
  - `max_iters`: The maximum number of iterations.
  - `boundaries`: Lower and upper boundaries for prevalences.
  - `boundaries_param`: Lower and upper boundaries for parameters.
  - `log`: Logical indicating whether to work with log weights.
  - `delete_induced_prior`: Logical indicating whether the induced prior density is to be deleted in the update of weights.
  - `mixture_samples`: The number of samples drawn from the weighted distribution to fit a new mixture to.
  - `df`: Degrees of freedom in the t-distributions, used yield a heavy tailed proposal.
  - `delta`, `sigma`, `breaks`: Options for density estimation used in the pseudo-likelihood and induced prior density.
- `seed`: Optional seed for the random number generator.
- `output_dir`: An optional string specifying the local directory where to save outputs after each iteration of the algorithm.
- `initial_amis_vals`: Optional list of intermittent outputs from a previous run.
  
  
### Defining a model function

The `amis` function expects its argument `model_func` to be a function with the 
following interface

```
observables <- model_func(seeds, params, n_tims)
```

- `params`: A matrix of parameter values (`double`) 
- `seeds`: A vector of seeds (`integer`)
- `n_tims`: Number of time points (`integer`)

Function `model_func` is expected to run the model for each pair
(`seed`, `params`) and return the corresponding values for the
observable (_e.g._ infection prevalence). `params` must be a matrix
with ncols equal to the dimension of the parameter space and the output
must be a matrix with ncols equal to the number of observed timepoints `n_tims`.

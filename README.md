# Description

This package provides an implementation of the Adaptive Multiple
Importance Sampling algorithm, as described in 

Integrating geostatistical maps and infectious disease transmission models using adaptive multiple importance sampling.
Renata Retkute, Panayiota Touloupou, María-Gloria Basáñez, T. Deirdre Hollingsworth and Simon E.F. Spencer (2021).
_Annals of Applied Statistics_, 15 (4), 1980-1998. DOI https://doi.org/10.1214/21-AOAS1486 

# Installation

Make sure you have the package [devtools](https://devtools.r-lib.org/)
installed. Then

```R
devtools::install_github("drsimonspencer/AMISforInfectiousDiseases")
```

# Usage

`amis` is the main function of the package. It takes a
geostatistical map, a transmission model and a prior distribution for parameters, 
and returns sampled parameters and their associated weights 
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
- `amis_params`: A list containing the parameters for the AMIS algorithm, such as:
  - `nsamples`: The number of sample parameters to draw at each iteration.
  - `boundaries`: Lower and upper boundaries for prevalences.
  - `mixture_samples`: The number of samples drawn from the weighted distribution to fit a new mixture to.
  - `target_ess`: The target effective sample size.
  - `max_iters`: The maximum number of iterations.
  - `delta`, `sigma`, `breaks`: Options for density estimation used in the Randon-Nikodym derivative.
  - `log`: logical indicating whether to work with log weights. 
  
## Defining a model function

The `amis` function expects its argument `model_func` to be a function with the 
following interface

```
observables <- model_func(seeds, params, n_tims)
```

- `params`: A matrix of parameter values (`double`) 
- `seeds`: A vector of seeds (`integer`)
- `n_tims`: Number of time points (`integer`)

Function `model_func` is expected to run the model for each pair
(`seed`, `parameter`) and return the corresponding values for the
observable (_e.g._ infection prevalence). `parameter` must be a matrix
with ncols equal to the dimension of the parameter space and the output
must be a matrix with ncols equal to the number of observed timepoints.

## Wrapping a model

If your model does not comply with the above interface, you can
provide a wrapper function. The following example illustrates how it
can be done in the case where the function implementing the model
reads and writes values in files and expects several additional
parameters to be specified:

```R
wrapped_model <- function(seeds, parameters) {
	# write input on disk
	input_file <- "beta_values.csv"
    write_model_input(seeds, parameters, input_file)
	
	# run model
    run_model(input_file, mda_file, output_file, infect_output,
               SaveOutput = F, OutSimFilePath = NULL,
               InSimFilePath = NULL)
			   
	# read results and return obsvervable
    res <- read.csv(output_file)
    return(100 * res[, dim(res)[2]])
  }
```

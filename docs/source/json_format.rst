JSON file format for model-balancing
====================================

.. _json:

Specifications
**************
The file format is based on JSON and has the following structure
(for a network with `Nr` reactions and `Nm` metabolites, and where `Ns` is the number of metabolic states
described in this dataset).

- `standard_concentration`: `str` - the standard concentration used for equilibrium constants
- `state_names`: `List[str]` - the names of the metabolic states (len = `Ns`)
- `network`: `dict`
    - `metabolite_names`: `List[str]` - the names of the metabolites (len = `Nm`)
    - `reaction_names`: `List[str]` - the names of the reactions (len = `Nr`)
    - `stoichiometric_matrix`: `array` - stoichiometric coefficients (shape = `(Nm, Nr)`)
    - `activation_matrix`: `array` - activation coefficients (shape = `(Nm, Nr)`)
    - `inhibition_matrix`: `array` - inhibition coefficients (shape = `(Nm, Nr)`)
- `kinetic_constants`: `dict`
    - `Keq`, `Kcatf`, `Kcatr`, `KM`, `KA`, or `KI`: `dict`
        - `unit`: `str` - the unit in which the values are given
        - `true`: `array` - true values in the case of artificial data (optional)
        - `prior_ln`: `dict`
            - `mean`: `array` - prior mean vector of log values
            - `cov`: `array` - prior covariance matrix of log values
        - `bounds_ln`: `dict`
            - `min`: `array` - lower bounds on log values
            - `max`: `array` - upper bounds on log values
        - `data_ln`: `dict`
            - `mean`: `array` - data mean vector of log values
            - `cov`: `array` - data covariance matrix of log values
        - `bounds`: `dict`
            - `min`: `array` - lower bounds
            - `max`: `array` - upper bounds
        - `combined`: `dict`
            - `geom_mean`: `array` - preposterior geometric mean vector
            - `mean_ln`: `array` - preposterior mean vector of log values
            - `cov_ln`: `array` - preposterior covariance matrix of log values
- `metabolite_concentrations` or `enzyme_concentrations`: `dict`
    - `unit`: `str` - the unit in which the values are given
    - `true`: `array` - true values in the case of artificial data (optional)
    - `prior_ln`: `dict`
        - `mean`: `array` - prior mean values for metabolite log-concentrations
        - `std`: `array` - prior std dev for metabolite log-concentrations
    - `bounds_ln`: `dict`
        - `min`: `array` - lower bounds on metabolite log-concentrations
        - `max`: `array` - upper bounds on metabolite log-concentrations
    - `data_ln`: `dict`
        - `mean`: `array` - data mean values for metabolite log-concentrations
        - `std`: `array` - data std dev for metabolite log-concentrations
    - `bounds`: `dict`
        - `min`: `array` - lower bounds
        - `max`: `array` - upper bounds
    - `combined`: `dict`
        - `geom_mean`: `array` - preposterior geom mean for metabolite concentrations
        - `geom_std`: `array` - preposterior geom std for metabolite concentrations
- `reaction_fluxes`: `dict`
    - `unit`: `str` - the unit in which the values are given
    - `true`: `array` - true values in the case of artificial data (optional)
    - `data`: `dict`
        - `mean`: `array` - flux data mean values
        - `std`: `array` - flux data std values

Examples
********

The directory (`examples/JSON`) stores examples for model balancing.

Note
****

- Activations and inhibitions are described by two separate (positive) matrices `activation_matrix` and 
  `inhibition_matrix`, where each one is of the shape `(Nm, Nr)`
- We follow the Matlab convention for ordering elements within a matrix. Therefore, the major index of the entries in 
  `kinetic_constants.KM`, `kinetic_constants.KA`, and `kinetic_constants.KI` (which are given in a sparse notation) is
  the metabolites and not the reactions (although rows correspond to the latter). In Python, the convention for ordering
  matrix elements is the opposite and therefore care must be taken.
- Fields called `combined` contain the preposterior distribution of the respective quantities


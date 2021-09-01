Prior distributions for model balancing
---------------------------------------

The prior table `cmb_prior.tsv` (in SBtab format) contains default information about priors, error bars, and constraints used in model balancing. The file is a modified version of the prior table `cmb_prior_ORIGINAL.tsv` used in parameter balancing (see [www.parameterbalancing.net] (www.parameterbalancing.net)).

## Change log

* replaced KM min 0.0000001 -> 0.00001
* replaced KM geom std 20 -> 5
* replaced Keq geom std 100 -> 1e8
* replaced concentration of enzyme geom mean 1 -> 0.001 (mM)
* replaced concentration geom mean 0.1 -> 0.3 (mM)
* replaced concentration min 10^-6 -> 10^-7 (mM)
* replaced product catalytic rate constant min 0.01 -> 0.005
* replaced substrate catalytic rate constant geom std 50 -> 1000
* replaced product catalytic rate constant geom std   50 -> 1000
* replaced DataGeomStd replaced for all quantities 1.2 -> 1.5

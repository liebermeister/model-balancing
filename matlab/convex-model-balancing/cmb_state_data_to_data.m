function data = cmb_state_data_to_data(state_data, cmb_options)

[nm,ns] = size(state_data.metabolite_data.Mean);

if isfield(state_data,'samples'),
  data.samples = state_data.samples;
else
  for it = 1:ns,
    data.samples{it,1} = ['S' num2str(it)];
  end
end

state_data.metabolite_data.Mean(state_data.metabolite_data.Mean < cmb_options.quantities.c.min) = cmb_options.quantities.c.min;
state_data.metabolite_data.Mean(state_data.metabolite_data.Mean > cmb_options.quantities.c.max) = cmb_options.quantities.c.max;

if length(state_data.metabolite_data.Std),
  [data.X.mean,data.X.std] = lognormal_normal_to_log(state_data.metabolite_data.Mean,state_data.metabolite_data.Std);
else
  data.X.mean = log(state_data.metabolite_data.Mean);
  data.X.std  = log(cmb_options.data_C_geom_std) * ones(nm,1);
end

data.E.mean = state_data.enzyme_data.Mean;
if length(state_data.enzyme_data.Std),
  data.E.std = state_data.enzyme_data.Std;
else
  data.E.std  = [cmb_options.data_E_geom_std-1] * data.E.mean;
end

data.V.mean = state_data.flux_data.Mean;
data.V.std  = state_data.flux_data.Std;
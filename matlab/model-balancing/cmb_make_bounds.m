function bounds = cmb_make_bounds(network,q_info,cmb_options, conc_min, conc_max)

eval(default('conc_min','[]','conc_max','[]'));
  
[nr,nm,nx,KM_ind,KA_ind,KI_ind,nKM,nKA,nKI] = network_numbers(network);

switch cmb_options.parameterisation
  case 'KV_KM_KA_KI',
    bounds.q_min = [log(cmb_options.quantities.KV.min) * ones(nr,1); ...
                    log(cmb_options.quantities.KM.min) * ones(nKM,1); ...
                    log(cmb_options.quantities.KA.min) * ones(nKA,1); ...
                    log(cmb_options.quantities.KI.min) * ones(nKI,1)];
    
    bounds.q_max = [log(cmb_options.quantities.KV.max) * ones(nr,1); ...
                    log(cmb_options.quantities.KM.max) * ones(nKM,1); ...
                    log(cmb_options.quantities.KA.max) * ones(nKA,1); ...
                    log(cmb_options.quantities.KI.max) * ones(nKI,1)];
  case 'Keq_KV_KM_KA_KI',
    nKeqind = length(q_info.q.index.Keq);
    bounds.q_min = [-log(cmb_options.quantities.Keq.max) * ones(nKeqind,1); ...
                    log(cmb_options.quantities.KV.min)   * ones(nr,1); ...
                    log(cmb_options.quantities.KM.min)   * ones(nKM,1); ...
                    log(cmb_options.quantities.KA.min)   * ones(nKA,1); ...
                    log(cmb_options.quantities.KI.min)   * ones(nKI,1)];
    
    bounds.q_max = [log(cmb_options.quantities.Keq.max) * ones(nKeqind,1); ...
                    log(cmb_options.quantities.KV.max)  * ones(nr,1); ...
                    log(cmb_options.quantities.KM.max)  * ones(nKM,1); ...
                    log(cmb_options.quantities.KA.max)  * ones(nKA,1); ...
                    log(cmb_options.quantities.KI.max)  * ones(nKI,1)];
  otherwise error('Ce_Options.Parameterisation not supported');
end

bounds.q_dep_min = [log(cmb_options.quantities.Kcatf.min) * ones(nr,1); ...
                    log(cmb_options.quantities.Kcatf.min) * ones(nr,1)];
bounds.q_dep_max = [log(cmb_options.quantities.Kcatf.max) * ones(nr,1); ...
                    log(cmb_options.quantities.Kcatf.max) * ones(nr,1)];

bounds.q_all_min = [bounds.q_min; ...
                    bounds.q_dep_min];
bounds.q_all_max = [bounds.q_max; ...
                    bounds.q_dep_max];

bounds.x_min = log(cmb_options.quantities.c.min * ones(nm,1));
bounds.x_max = log(cmb_options.quantities.c.max * ones(nm,1));

if length(conc_min),
  bounds.x_min = log(conc_min);
  bounds.x_max = log(conc_max);
end

bounds.v_min = -cmb_options.quantities.v.max * ones(nr,1);
bounds.v_max =  cmb_options.quantities.v.max * ones(nr,1);

bounds.e_min =  zeros(nr,1);
bounds.e_max =  cmb_options.quantities.e.max * ones(nr,1);

bounds.a_forward_min = cmb_options.quantities.Aforward.min * ones(nr,1);
bounds.a_forward_max = cmb_options.quantities.Aforward.max * ones(nr,1);

bounds.conc_min = conc_min;
bounds.conc_max = conc_max;

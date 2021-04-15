function [V, Vstd] = cmb_project_fluxes(V,Vstd,network,flag_project_fluxes)

% [V, Vstd] = cmb_project_fluxes(V,Vstd,network,flag_project_fluxes)
%
% V, Vstd: Flux mean and std (matrices)
% flag_project_fluxes: 'all', 'missing', 'none'

eval(default('flag_project_fluxes','''all'''));


switch flag_project_fluxes,
  case {'all','missing'},
    display(sprintf('Replacing %s fluxes by projections',flag_project_fluxes));
    %% Complete flux data by using projection (MOVE THIS ELSEWHERE???)
    [nr,ns]  = size(V);
    my_V     = V;
    my_V_std = Vstd;
    v_sign   = sign(mean(sign(my_V),2)+0.1);
    for it = 1:ns,
      [my_V_proj(:,it), ~, my_V_proj_std(:,it)] = project_fluxes(network.N, find(network.external), my_V(:,it), my_V_std(:,it),v_sign);
    end
end

switch flag_project_fluxes,
  case 'all',
    V = my_V_proj;
    Vstd  = my_V_proj_std;
  case 'missing',
    ind = find(~isfinite(V));
    V(ind) = my_V_proj(ind);
    Vstd(ind)  = my_V_proj_std(ind);
  end
end

function [ estCap, varargout ] = estimate_capacitance_Outer_fn( tslist_cell, Vlist_array, ...
    gelist, gilist, capoutparam)
%  Return a single number
%       

% Default parameters, which you can overwrite in "capoutparam"

% Values for testing capacitance and averaging. 
% Why are these different? You may want to make sure averaging region is properly chosen 
Vrest_min = -80;
Vrest_max = -50;
Vavg_min = -67;
Vavg_max = -62;
tau_ref  = 0.1;  % In seconds

if (isfield(capoutparam,'Vrest_min')); Vrest_min = capoutparam.Vrest_min; end;
if (isfield(capoutparam,'Vrest_max')); Vrest_max = capoutparam.Vrest_max; end;
if (isfield(capoutparam,'Vavg_min')); Vavg_min = capoutparam.Vavg_min; end;
if (isfield(capoutparam,'Vavg_max')); Vavg_max = capoutparam.Vavg_max; end;
if (isfield(capoutparam,'tau_ref')); tau_ref = capoutparam.tau_ref; end;

% Now: convert voltage to V
Vrest_array     = Vrest_min:Vrest_max;        % in mV
Vrest_array     = Vrest_array/1000;   % in V
Vavg_min    = Vavg_min/1000;
Vavg_max    = Vavg_max/1000;

% Parameter structure for INNER function
capparam = [];

% Parameters that need to be copied over to inner structure
if (isfield(capoutparam,'dt')); capparam.dt = capoutparam.dt; end;
if (isfield(capoutparam,'Ve')); capparam.Ve = capoutparam.Ve; end;
if (isfield(capoutparam,'Vi')); capparam.Vi = capoutparam.Vi; end;
if (isfield(capoutparam,'Aexc')); capparam.Aexc = capoutparam.Aexc; end;
if (isfield(capoutparam,'Ainh')); capparam.Ainh = capoutparam.Ainh; end;

% EXCLUDE some reasonable time after spike
capparam.tau_ref = tau_ref;  %in s

Var_array_all = [];
estCap_all = []; Vrest_keep =[];

for j1=1:length(Vrest_array)
   capparam.Vrest = Vrest_array(j1);
   [estCap,blah]=estimate_capacitance_fn( tslist_cell, Vlist_array, gelist, gilist, capparam ); 
   
   if (~isempty(estCap))
       Vrest_keep = [Vrest_keep Vrest_array(j1)];
       estCap_all = [estCap_all estCap];
       Var_array_all = [Var_array_all blah(:,2)];
   end
end

if (nargout > 1)
   % Send back more detailed information...
   moreInfo.Vrest_keep = Vrest_keep;
   moreInfo.estCap_all = estCap_all;
   moreInfo.Var_array_all = Var_array_all;
   varargout{1} = moreInfo;
end

%
% Badel et al says: look near reversal potential.
%       (However, our estimate starts to fluctuate near that level; makes
%       me leary!)
to_avg_ind = find(Vrest_keep >= Vavg_min & Vrest_keep <= Vavg_max);

estCap = mean(estCap_all(to_avg_ind));     % in nF



end


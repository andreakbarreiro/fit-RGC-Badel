function [params,varargout] = fit_f_fminsearch(V, f, varargin)

init_cond = [1/3 -88.5 4 -31.5];

% Optional: pass in IC
if (nargin > 2)
    init_cond = varargin{1};
end

options=[];   
%options=optimset('Display','iter');    % While debugging

options = optimset(options,'MaxIter',1000,'MaxFunEvals',2000);

[params,fval,fflag,outpt] = fminsearch(@f_fit, init_cond, options, V, f);

if (nargout > 1)
    varargout{1} = fflag;
end
if (nargout > 2)
    varargout{2} = outpt;
end
if (nargout > 3)
    varargout{3} = fval;
end

approx = params(1)*(params(2)-V+params(3)*exp((V-params(4))/params(3)));
%figure;plot(V, f); hold on;plot(V,approx,'r');
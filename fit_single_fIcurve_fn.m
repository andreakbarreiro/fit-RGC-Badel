function [ tscore, ffit_many_sendback, ffit_flag_sendback ] = fit_single_fIcurve_fn( Varray, farray )
% fit single fIcurve
%   
%    Varray:        voltage values (in mV)
%    farray:        f values
% 

% How many repeats?
nTry = 10;

ffit_many      = zeros(4,nTry);
ffit_good_many  = zeros(1,nTry);
ffit_score_many = zeros(1,nTry);

% Background
init_bg = [1/3 -88.5 4 -31.5];

for j1=1:nTry
    initg = init_bg .* (ones(1,4) + 0.2*randn(1,4));
   
    
    [ptemp,fflag] = fit_f_fminsearch(Varray,farray,initg);
    
    ffit_many(:,j1)         = ptemp;
    ffit_good_many(j1)      = fflag;
    ffit_score_many(j1)     = f_fit(ptemp,Varray,farray);

end

% Order in ascending order...
[tscore,ScInd]=sort(ffit_score_many,'ascend');

ffit_many_sendback = zeros(size(ffit_many));
for k2=1:nTry
   ffit_many_sendback(:,k2) = ffit_many(:,ScInd(k2)); 
end

ffit_flag_sendback = ffit_good_many(ScInd);


end


function [ tslist,tplist,peaklist ] = get_spike_times_fn( Vlist, spikeparam )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % Default parameters
    dt    = .0001;
    % Official "threshold"
    Vth_ck  = 0;
    % This is where we record the crossing
    Vth_cr  = -0.02; 
    % Refractory period
    tau_ref = 0.002;

    % Process any parameters
    if (isfield(spikeparam,'dt')); dt = spikeparam.dt; end;
    if (isfield(spikeparam,'Vth_ck')); Vth_ck = spikeparam.Vth_ck; end;
    if (isfield(spikeparam,'Vth_cr')); Vth_cr = spikeparam.Vth_cr; end;
    if (isfield(spikeparam,'tau_ref')); tau_ref = spikeparam.tau_ref; end;

    nTRef = round(tau_ref/dt);

    totalT  = (length(Vlist)-1)*dt;
    tlist = 0:dt:totalT;

    tslist = [];
    tplist = [];
    peaklist = [];
    n = 1;
    while (n < numel(tlist))
        if (Vlist(n) > Vth_ck)    % a spike has occurred
            tslist = [tslist tlist(n)];
    
            % Find peak
            nflag=0;j=n;
            while (nflag==0)
                if (Vlist(j+1)<Vlist(j))   % Assume peak occurs right before voltage 
                                   % decreases
                    tplist = [tplist tlist(j)];
                    peaklist = [peaklist Vlist(j)];
                    nflag = 1;
                else
                    j=j+1;
                end
            end

            % Now, advance index to reflect absolute refraction period;
            %      don't want to pick up "downstroke" of the spike
            n = n+nTRef;  % # of time steps 
        else
            n=n+1;
        end
    end

end % function


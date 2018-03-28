function [Var_array] = collectData( tslist_cell, Vlist_array, gelist, gilist, capparam )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % Default parameters
    Vrest = -70/1000;           % Resting voltage
    dt    = .0001;              % Time step (in s)


    Ve = 0;                     % E reversal potential
    Vi = -.08;                  % I reversal potential
    Aexc = 30;                  % Normalization factor for Aexc
    Ainh = 40;

    tau_ref = 0.003;            % Time to discard after spike

    % Process any parameters
    if (isfield(capparam,'dt')); dt = capparam.dt; end;
    if (isfield(capparam,'Vrest')); Vrest = capparam.Vrest; end;
    if (isfield(capparam,'Ve')); Ve = capparam.Ve; end;
    if (isfield(capparam,'Vi')); Vi = capparam.Vi; end;
    if (isfield(capparam,'Aexc')); Aexc = capparam.Aexc; end;
    if (isfield(capparam,'Ainh')); Ainh = capparam.Ainh; end;
    if (isfield(capparam,'tau_ref')); tau_ref = capparam.tau_ref; end;

    % Did I actually pass a cell of spike times?
    tlistcell_flag = 0;
    if (iscell(tslist_cell)); tlistcell_flag = 1; end;


    % Pick voltage at which to measure capacitance (in Volts)
    Vmin  = Vrest-(.5/1000);
    Vmax  = Vrest+(.5/1000);

    % Potential Capacitances; in nF (Fred estimates 50 pF)
    Cap_array = [1:100 110:10:300]/1000;
    Var_array = zeros(size(Cap_array));

    totalT  = (length(Vlist_array)-1)*dt;
    tlist = 0:dt:totalT;

    % For each voltage data stream, find spike times and exclude voltage
    % measurements right after spike.
    nTRef = round(tau_ref/dt);

    dVdtlist = [];
    Iapplist = [];

    % How many spike time lists to process?
    if (tlistcell_flag)
        nList = length(tslist_cell);
    else
        nList=1;
    end

    for k=1:nList
    %
        if (tlistcell_flag)
            tslist = tslist_cell{k};
        else
            tslist = tslist_cell;
        end
        Vlist  = Vlist_array(k,:);

        win_labels = zeros(1,length(tlist)-1);

        % Exclude "tau_ref" ms after spike time
        tpad=2*dt;   % To exclude BEFORE next spike
        for j=1:length(tslist)-1
            ts=tslist(j);           % Current spike time
            tsnext=tslist(j+1);     % Next spike time

            win_labels(find( ts < tlist & tlist < min(ts+tau_ref,tsnext-tpad) ))=1;
        end

        % Needed data 
        dVdt_temp      = diff(Vlist)/dt; Vlist = Vlist(1:end-1);
        Iapplist_temp  = Aexc*gelist.*(Ve-Vlist)+Ainh*gilist.*(Vi-Vlist);
        %Iapplist_temp  = Iapplist_temp(1:end-1);

        %tempind = find(Vlist >= Vmin & Vlist <= Vmax); 

        okind     = find(Vlist >= Vmin & Vlist <= Vmax & win_labels==0);

        Iapplist  = [Iapplist Iapplist_temp(okind)];
        dVdtlist  = [dVdtlist dVdt_temp(okind)];
    end

    for i=1:numel(Cap_array)
       C            = Cap_array(i);
       data         = (Iapplist - C*dVdtlist)/C;
       %data         = (Iapplist/C) - dVdtlist;
       Var_array(i) = var(data);
    end
end


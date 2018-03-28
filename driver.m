% ATTENTION: CALL AFTER get_spike_times_driver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Main Driver                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IF you want to override default parameters
capoutparam = [];

buildPlot(tslist_cell, Vlist_array, gelist, gilist, capoutparam);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Function Declaration                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = buildPlot( tslist_cell, Vlist_array, ...
    gelist, gilist, capoutparam )
%buildPlot This function collects data from the spike times and plots it

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             Boilerplate Code                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Default parameters, which you can overwrite in "capoutparam"

    % Values for testing capacitance and averaging. 
    % Why are these different? You may want to make sure averaging region is 
    % properly chosen 
    Vrest_min = -80;
    Vrest_max = -50;
    Vavg_min = -67;
    Vavg_max = -62;
    tau_ref  = 0.1;  % In seconds

    if (isfield(capoutparam,'Vrest_min'))
        Vrest_min = capoutparam.Vrest_min; 
    end;
    if (isfield(capoutparam,'Vrest_max')) 
        Vrest_max = capoutparam.Vrest_max; 
    end;
    if (isfield(capoutparam,'Vavg_min')) 
        Vavg_min = capoutparam.Vavg_min; 
    end;
    if (isfield(capoutparam,'Vavg_max')) 
        Vavg_max = capoutparam.Vavg_max; 
    end;
    if (isfield(capoutparam,'tau_ref'))
        tau_ref = capoutparam.tau_ref; 
    end;

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
    
    X = Vrest_array;
    Y = [1:100 110:10:300]/1000;
    Z = zeros(31,120);
    tempZ = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                               Collect Data                               %
    % This section calls the collectDate() function to get the Variance of the %
    % data in order to plot it later                                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:length(Vrest_array)
       capparam.Vrest = Vrest_array(i);
       tempZ =collectData(tslist_cell, Vlist_array, gelist, gilist, capparam); 
       
       Z(i, :) = tempZ;
    end
    
    %---Checking the sizes of everything---%
    size(X)
    size(Y)
    size(Z)
    %--------------------------------------%
    
    % <WHAT DOES THIS DO? and why 3?>
    view(3)
    
    % Create a grid of points to plot because plot3 needs equal sized matrices
    [Xmat,Ymat]=meshgrid(X,Y);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    Creating the plot of the minimums                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    B = ones(120, 1);
    B = B.*0.02;
    [Amat, Bmat] = meshgrid(X,B);
    C = Z';
    
    % Finding the minimum of C
    [M,I] = min(Z');
    [ZZ, ZY] = meshgrid(M,M);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Plotting everything and Labeling it                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    plot3(Xmat,Ymat,Z,'*');
    hold on;
    plot3(X, Ymat(I), M, 'ro');
    
    %Labels
    xlabel('V');
    ylabel('C_{est}');
    zlabel('Variance');
    set(gca,'FontSize',16);
    zlim([0,5])
end  % function buildPlot()
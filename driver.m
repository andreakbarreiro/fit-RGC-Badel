%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Main Driver                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the workspace is empty, run get_spike_times_driver
if ( ~exist('dfile', 'var') )
    get_spike_times_driver
end

% Override default parameters in organizeData()
organizeParam = [];

% Organize and collect data
[success, dataArrays] = organizeData(tslist_cell, Vlist_array,  ...
                                     gelist, gilist, organizeParam);

% Unpack dataArrays
Vrest_array = dataArrays.Vrest;
Cap_array = dataArrays.Cap;
Var_array = dataArrays.Var;
Mean_array = dataArrays.Mean;
numV_array = dataArrays.numV;

% Finding the minimum of the variances
[MinVar,IndexVar] = min(Var_array');

% Store the value of the capacitance where the minimum is
minCapArray = Cap_array(IndexVar);

% Build Plots of data 
buildPlot(Vrest_array, Cap_array, Var_array, Mean_array, MinVar, IndexVar)

% Calculate the intercept
meanInter = findIntercept(Vrest_array, Mean_array)

% Calculate an estimated true capacitance
estTrueCap = dataToGaussian(Vrest_array, minCapArray, meanInter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Function Declarations                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [runFlag, varargout] = organizeData( tslist_cell, Vlist_array, ...
                            gelist, gilist, ODParam )
% organizeData Returns arrays of the variance and mean of the data
%   Input
%       tslist_cell = spike times from get_spike_times_driver
%       Vlist_array = Voltages
%       gelist      = Conductances for variance calculation
%       gilist      = Conductances for variance calculation
%       ODParam  = Array for parameters to run with
%
%   Output
%       runFlag     = A zero on a successful function execution, one otherwise
%       varargout   = Vrest, Capacitance, Variance, Mean, and numV arrays 
%                     in that order

    % Default parameters
    Vrest_min = -80;
    Vrest_max = -50;
    Vavg_min = -67;
    Vavg_max = -62;
    tau_ref = 0.1;  % In seconds

    % Checking if ODParam overrides default parameters
    if (isfield(ODParam,'Vrest_min'))
        Vrest_min = ODParam.Vrest_min; 
    end;
    if (isfield(ODParam,'Vrest_max')) 
        Vrest_max = ODParam.Vrest_max; 
    end;
    if (isfield(ODParam,'Vavg_min')) 
        Vavg_min = ODParam.Vavg_min; 
    end;
    if (isfield(ODParam,'Vavg_max')) 
        Vavg_max = ODParam.Vavg_max; 
    end;
    if (isfield(ODParam,'tau_ref'))
        tau_ref = ODParam.tau_ref; 
    end;
    
    % Initialize Vrest_array
    Vrest_array = Vrest_min:Vrest_max;
    
    % Convert voltage to Volts
    Vrest_array = Vrest_array/1000;
    Vavg_min = Vavg_min/1000;
    Vavg_max = Vavg_max/1000;

    % Initialize collectData() parameter array
    CDParam = [];

    % Parameters that need to be copied over to inner structure
    if (isfield(ODParam,'dt')); CDParam.dt = ODParam.dt; end;
    if (isfield(ODParam,'Ve')); CDParam.Ve = ODParam.Ve; end;
    if (isfield(ODParam,'Vi')); CDParam.Vi = ODParam.Vi; end;
    if (isfield(ODParam,'Aexc')); CDParam.Aexc = ODParam.Aexc; end;
    if (isfield(ODParam,'Ainh')); CDParam.Ainh = ODParam.Ainh; end;

    % Exclude some reasonable time after spike
    CDParam.tau_ref = tau_ref;  %in s

    % Initialize storage arrays
    Var_array_all = [];
    estCap_all = []; 
    Vrest_keep =[];
    varData = [];
    meanData = [];
    
    % Set up the Capacitance array
    Cap_array = [1:100 110:10:300]/1000;
    
    % Set Z axis
    Var_array = zeros(31,120);
    Mean_array = zeros(31,120);
    
    numV_array = zeros(size(Vrest_array));

    % Call collectData() collect the desired data
    for i=1:length(Vrest_array)
        % Grab a slice of Vrest_array
        CDParam.Vrest = Vrest_array(i);
        
        % Collect data 
        [varData, meanData, numV] = collectData(tslist_cell, Vlist_array, ...
                                               gelist, gilist, CDParam); 
        % Append new data
        Var_array(i, :) = varData;
        Mean_array(i, :) = meanData;
        
        % Unused
        numV_array(i) = numV;
    end
    
    % Send back more detailed info
    if (nargout > 1)
        dataArrays.Vrest = Vrest_array;
        dataArrays.Cap = Cap_array;
        dataArrays.Var = Var_array;
        dataArrays.Mean = Mean_array;
        dataArrays.numV = numV_array;
        varargout{1} = dataArrays;
    end
    
    runFlag = 0;
    
end %function organizeData()

function [Var_array, Mean_array, numV] = collectData( tslist_cell, Vlist_array, gelist, gilist, capparam )
%collectData A helper function to buildPlot() that runs through data and pulls
%            out relevant info
%   Input
%       tslist_cell = spike times from get_spike_times_driver
%       Vlist_array = Voltages
%       gelist      = Conductances for variance calculation
%       gilist      = Conductances for variance calculation
%       ODCapParam  = Array for parameters to run with
%
%   Output
%       Var_array    = Array of variances
%       Mean_array   = Array of means
%       numV         = Length of dVdtlist 

    % Default parameters
    Vrest = -70/1000;      % Resting voltage
    dt = .0001;            % Time step (in s)
    Ve = 0;                % E reversal potential
    Vi = -.08;             % I reversal potential
    Aexc = 30;             % Normalization factor for Aexc
    Ainh = 40;
    tau_ref = 0.003;       % Time to discard after spike
    
    % Initialize variables
    tlistcell_flag = 0;    % Flag for cell of spike times
    
    % Process any parameters
    if (isfield(capparam,'dt')); dt = capparam.dt; end;
    if (isfield(capparam,'Vrest')); Vrest = capparam.Vrest; end;
    if (isfield(capparam,'Ve')); Ve = capparam.Ve; end;
    if (isfield(capparam,'Vi')); Vi = capparam.Vi; end;
    if (isfield(capparam,'Aexc')); Aexc = capparam.Aexc; end;
    if (isfield(capparam,'Ainh')); Ainh = capparam.Ainh; end;
    if (isfield(capparam,'tau_ref')); tau_ref = capparam.tau_ref; end;

    % Checking a cell of spike times was passed in as a parameter
    if (iscell(tslist_cell)); tlistcell_flag = 1; end;

    % Set range of voltages (in Volts) to measure capacitance
    Vmin  = Vrest-(.5/1000);
    Vmax  = Vrest+(.5/1000);

    % Potential Capacitances; in nF (Fred estimates 50 pF)
    % Initialize storage arrays
    Cap_array = [1:100 110:10:300]/1000;
    Var_array = zeros(size(Cap_array));
    Mean_array = zeros(size(Cap_array));
    dVdtlist = [];
    Iapplist = [];
    
    % ???
    % Set total time interval, discretize time
    totalT  = (length(Vlist_array)-1)*dt;
    tlist = 0:dt:totalT;

    % For each voltage data stream, find spike times and exclude voltage
    % measurements right after spike.
    nTRef = round(tau_ref/dt);

    % If we have a cell of spike times, get its length
    if (tlistcell_flag)
        nList = length(tslist_cell);
    else
        nList=1;
    end

    % For every item in the cell of spike times
    for k=1:nList
        
        if (tlistcell_flag)
            tslist = tslist_cell{k};
        else
            tslist = tslist_cell;
        end
        
        % Grab a slice of Vlist_array
        Vlist  = Vlist_array(k,:);

        % Initialize storage for __peak times__
        win_labels = zeros(1,length(tlist)-1);

        % Exclude "tau_ref" ms after spike time and before next spike
        tpad=2*dt;
        
        for j=1:length(tslist)-1
            % Current spike time
            ts=tslist(j);
            
            % Next spike time
            tsnext=tslist(j+1);

            % Mark a window of time after a spike to exclude from calculation
            win_labels( ...
                find( ts < tlist  ...
                      & tlist < min(ts+tau_ref,tsnext-tpad)  ...
                     ) ...
            ) = 1;
        end

        % Calculate needed data 
        dVdt_temp = diff(Vlist)/dt; Vlist = Vlist(1:end-1);
        Iapplist_temp = Aexc*gelist.*(Ve-Vlist)+Ainh*gilist.*(Vi-Vlist);

        % 
        okind = find(Vlist >= Vmin & Vlist <= Vmax & win_labels==0);
        
        % Append current iteration of data
        Iapplist  = [Iapplist Iapplist_temp(okind)];
        dVdtlist  = [dVdtlist dVdt_temp(okind)];
    end
    
    % Return the length of dVdtlist 
    numV = length(dVdtlist);
    
    for i=1:numel(Cap_array)
        % Intermediate calculations
        C            = Cap_array(i);
        data         = (Iapplist - C*dVdtlist)/C;
        
        % Storing variance and mean of the data
        Var_array(i) = var(data);
        Mean_array(i) = mean(data);
    end
end % end function collectData()

function [] = buildPlot(Vrest_array, Cap_array, Var_array, ...
                        Mean_array, MinVar, IndexVar)
% buildPlot This function collects data from the spike times and plots it
%   Input
%       Vrest_array = Array of resting voltages 
%       Cap_array   = Array of capacitances 
%       Var_array   = Array of variances 
%       Mean_array  = Array of means
%       MinVar      = Array of the values of minimized capacitances 
%       IndexVar    = Array of the indexes of minimized capacitances
%
%   Output
%       3 plots

    % Create a grid of points to plot because plot3 needs equal sized matrices
    [Xgrid,Ygrid] = meshgrid(Vrest_array,Cap_array);
    
    % Figure 1: Voltage and Capacitance versus Variance of the data
    plot3(Xgrid,Ygrid,Var_array,'*');
    
    % Allow plots to be added to
    hold on;
    
    % Add to Figure 1 the minimum of each Capacitance-Variance curve
    plot3(Xgrid, Ygrid(IndexVar), MinVar, 'ro');
    
    % Labels for Figure 1
    title('Variance Plot');
    xlabel('V');
    ylabel('C_{est}');
    zlabel('Variance');
    set(gca,'FontSize',16);
    zlim([0,5])
    
   % Figure 2: Voltage and Capacitance versus Mean of the data
    figure;
    plot3(Ygrid',Xgrid',Mean_array','*');
    
    % Labels for Figure 2
    title('Mean Plot');
    xlabel('C_{est}');
    ylabel('V');
    zlabel('Mean');
    set(gca,'FontSize',16);
    
    % Figure 3: 
    figure;
    plot(1000*Vrest_array', Mean_array(:,1), '*')
    hold on;
    
    % Plotting the middle set of voltages
    plot(1000*Vrest_array(14:end-8)', 0, '+')
    
    % Labels for Figure 3
    title('Finding the intercept at a fixed voltage')
    xlabel('Voltage')
    ylabel('Mean')
    set(gca,'FontSize',16);
    
end % function buildPlot()

function [meanInter] = findIntercept(Vrest_array, Mean_array)
% buildPlot This function collects data from the spike times and plots it
%   Input
%       Vrest_array = Voltages
%       Mean_array  = Mean array from organizeData()
%
%   Output
%       meanInter = Mean of the calculated intercepts

    % Isolate data from fixed voltages
    fixedVolt1 = Mean_array(:,1);
    fixedVolt2 = Mean_array(:,2);
    fixedVolt3 = Mean_array(:,3);

    % Linear fits to points near the intercept
    P = polyfit(1000*Vrest_array(10:end-9)',fixedVolt1(10:end-9),1);
    Q = polyfit(1000*Vrest_array(10:end-9)',fixedVolt2(10:end-9),1);
    R = polyfit(1000*Vrest_array(10:end-9)',fixedVolt3(10:end-9),1);

    % Calculate individual intercepts
    pval = -P(2) / P(1);
    qval = -Q(2) / Q(1);
    rval = -R(2) / R(1);

    % Package together the intercepts
    inter = [pval, qval, rval];

    % return the mean of the intercepts
    meanInter = mean(inter);

end % function findIntercept()

function [estTrueCap] = dataToGaussian(Vrest_array, minCap, meanInter)
% dataToGaussian Estimates capacitance based on normal distribution of points
%   Input
%       Vrest_array = Array of Voltage Data
%       minCap      = Array of minimum capacitance values
%       meanInter   = Mean from findIntercept()
%
%   Output
%       estTrueCap  = Estimated True Capacitance based on minimized variance

    % Find the normal distribution with mean=meanInter2 and SD=1
    nodeWeights = normpdf(1000*Vrest_array(14:end-8), meanInter, 1);

    % Estimate capacitance from the middle 10 voltages
    estTrueCap = dot(minCap(14:end-8), nodeWeights)
    
end % function dataToGaussian()